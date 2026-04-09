import subprocess
import argparse
from sqlalchemy import create_engine
import pandas as pd
from datetime import datetime, timezone
from uuid import uuid4
from zoneinfo import ZoneInfo

parser = argparse.ArgumentParser(
        description="Upload demultiplex stats, run info, and sample sheet metadata to SMDB."
    )
parser.add_argument("-c", "--path_to_demultiplex_stats", required=True, help="Path to demultiplex stats CSV")
parser.add_argument("-r", "--path_to_run_info", required=True, help="Path to run info XML")
parser.add_argument("-x", "--path_to_sample_sheet", required=True, help="Path to sample sheet CSV")
parser.add_argument("-n", "--db_name", required=True, help="Database name")
parser.add_argument("-s", "--schema_name", required=True, help="Target schema name")
parser.add_argument("-u", "--db_user", required=True, help="Database user")
parser.add_argument("-p", "--db_password", required=True, help="Database password")
parser.add_argument("-d", "--db_host", required=True, help="Database host")
parser.add_argument("-o", "--db_port", required=True, type=int, help="Database port")
parser.add_argument("-t", "--table_name", required=True, help="Target table name")
parser.add_argument("-e", "--send_upload_receipts_to", required=True, help="Space-separated emails for upload receipts")

args = parser.parse_args()

def upload_demultiplex_stats(path_to_demultiplex_stats, path_to_run_info, path_to_sample_sheet, db_name, schema_name, db_user, db_password, db_host, db_port, table_name, send_upload_receipts_to):
    """
    Connects to a PostgreSQL database, queries a view, and exports the results to a TSV file.

    Args:
        path_to_demultiplex_stats (str): Path to the demultiplex stats CSV file.
        path_to_run_info (str): Path to the run info CSV file.
        path_to_sample_sheet (str): Path to the sample sheet CSV file.
        db_name (str): Database name.
        db_user (str): Database user.
        db_password (str): Database password.
        db_port (int): Database port (default is 5432).
        table_name (str): Name of the PostgreSQL table to export.
        output_file (str): Path to the output TSV file.
    """
   
    dmux_stats = pd.read_csv(path_to_demultiplex_stats)
    run_info = pd.read_xml(path_to_run_info, parser='etree')
    sample_sheet = pd.read_csv(path_to_sample_sheet)
    
    sample_sheet = (
        sample_sheet.iloc[0:17]
        .dropna(how="all", axis=1)
        .dropna(how="any", axis=0)
        .rename(columns={"[Header]": "Attribute", "Unnamed: 1": "Value"})
        .set_index("Attribute")
    )
    
    send_upload_receipts_to =  send_upload_receipts_to.replace(" ", "; ")
    upload_sheets = str(path_to_demultiplex_stats) + "; " + str(path_to_run_info) + "; " + str(path_to_sample_sheet)
    timestamp = datetime.now(ZoneInfo("Europe/Copenhagen")) 
    seq_run_id = run_info.at[0, 'Id']
    seq_run_number = run_info.at[0, 'Number']
    seq_machine_id = run_info.at[0, 'Instrument']
    flowcell_position = run_info.at[0, 'Id'].split('_')[3][0]
    flowcell_id = run_info.at[0, 'Flowcell']
    
    unformatted_seq_date = run_info.at[0, 'Date']
    dt = datetime.strptime(unformatted_seq_date, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    formatted_date = dt.strftime("%Y-%m-%d")
    
    dmux_stats['database_insert_by'] = send_upload_receipts_to
    dmux_stats['upload_sheet'] = upload_sheets
    dmux_stats['database_insert_datetime_utc'] = timestamp
    dmux_stats['upload_uuid'] = uuid4()
    dmux_stats['sequencing_run_id'] = seq_run_id
    dmux_stats['sequencing_run_number'] = seq_run_number
    dmux_stats['sequencing_machine_id'] = seq_machine_id
    dmux_stats['sequencing_run_date'] = str(formatted_date)
    dmux_stats['sequencing_tube_tag'] = None
    dmux_stats['flowcell_position'] = flowcell_position    
    dmux_stats['flowcell'] = flowcell_id
    dmux_stats.loc[dmux_stats["SampleID"] == "Undetermined", "Index"] = "unknown"
    
    # TODO: Uncomment
    num_lanes = len(dmux_stats['Lane'].unique())
    assert 0 < num_lanes < 9, f"Expected 1-8 lanes in the demultiplex stats file, but got {num_lanes}"
    
    # Add sequencing pool based on lane
    pool_lanes = [sample_sheet.loc['PoolLane1', 'Value'],
                  sample_sheet.loc['PoolLane2', 'Value'],
                  sample_sheet.loc['PoolLane3', 'Value'],
                  sample_sheet.loc['PoolLane4', 'Value'],
                  sample_sheet.loc['PoolLane5', 'Value'],
                  sample_sheet.loc['PoolLane6', 'Value'],
                  sample_sheet.loc['PoolLane7', 'Value'],
                  sample_sheet.loc['PoolLane8', 'Value']
                  ]
    
    dmux_stats['sequencing_tube_tag'] = dmux_stats['Lane'].apply(lambda x: pool_lanes[x - 1])
    
    ENGINE = create_engine(f"postgresql://{db_user}:{db_password}@{db_host}:{db_port}/{db_name}")    
    
    q = 'select column_name_db, column_name_sheet from name_maps.column_names where table_id = 9;'

    renamer = pd.read_sql(q, ENGINE)
    rename_dict = dict(zip(renamer['column_name_sheet'], renamer['column_name_db']))
    dmux_stats = dmux_stats.rename(columns=rename_dict)
    
    dmux_stats.to_sql(table_name, ENGINE, schema=schema_name, if_exists="append", index=False)

    email_cmd = [
        "mail",
        "-s",
        "Sequencing stats successfully uploaded to SMDB",
        "-a",
        str(path_to_demultiplex_stats),
        "-a",
        str(path_to_run_info),
        "-a",
        str(path_to_sample_sheet)
    ]
    email_cmd.extend(send_upload_receipts_to.split("; "))
    subprocess.run(
        email_cmd,
        input="The appended sequencing stats has been successfully uploaded to SMDB.",
        text=True,
        check=False,
    )
    

upload_demultiplex_stats(
    path_to_demultiplex_stats=args.path_to_demultiplex_stats,
    path_to_run_info=args.path_to_run_info,
    path_to_sample_sheet=args.path_to_sample_sheet,
    db_name=args.db_name,
    schema_name=args.schema_name,
    db_user=args.db_user,
    db_password=args.db_password,
    db_port=args.db_port,
    table_name=args.table_name,
    db_host=args.db_host,
    send_upload_receipts_to=args.send_upload_receipts_to
    )
    

    
