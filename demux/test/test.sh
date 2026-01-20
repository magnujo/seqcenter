
python3 ../cross_contamination.py --loglevel ERROR --index-counts 20251013_A00706_0986_AHHTKKDSXF.Index_Hopping_Counts.csv --index-known ../eDNA_index_list_UDP097-UDP288_UDI001-UDI096_250807.txt --lanes 1 --rpm-warn 100 --out-prefix /tmp/cc_test1
python3 ../cross_contamination.py --loglevel ERROR --index-counts 20251013_A00706_0986_AHHTKKDSXF.Index_Hopping_Counts.csv --index-known ../eDNA_index_list_UDP097-UDP288_UDI001-UDI096_250807.txt --lanes 2 --rpm-warn 100 --out-prefix /tmp/cc_test2
python3 ../cross_contamination.py --loglevel ERROR --index-counts 20251013_A00706_0986_AHHTKKDSXF.Index_Hopping_Counts.csv --index-known ../eDNA_index_list_UDP097-UDP288_UDI001-UDI096_250807.txt --lanes 3,4 --rpm-warn 100 --out-prefix /tmp/cc_test34
md5sum -c test.md5
rm /tmp/cc_test*
