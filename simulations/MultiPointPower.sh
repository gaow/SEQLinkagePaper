wsqlite MultiPoint.sqlite3 -s CHP
wsqlite MultiPoint.sqlite3 -s MultiPoint
wsqlite MultiPoint.sqlite3 "select * from multipoint where fn like '%FAM20%gene0%'"
wsqlite MultiPoint.sqlite3 "select fam_size,prop1,plod1 from CHP where fam_size=20"
# wsqlite PowerCalc.sqlite3 "select fam_size,prop1,plod1 from CHP where fam_size=21 and moi='compound_recessive' and ahet=1"
wsqlite PowerCalc.sqlite3 "select fam_size,prop1,plod1 from SNV where fam_size=21 and moi='compound_recessive' and ahet=1"
