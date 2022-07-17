all:
	cd doxygen && make
	cd docbook && make

WWW=/home/pae/praga/WWW/pangolin/doc/
FOLDER_DOXY=pangolin_ref
FOLDER_DOCB=user_manual

# Upload to cerfacs website
upload:
	mkdir -p ${WWW}/${FOLDER_DOXY}
	rsync -avz doxygen/${FOLDER_DOXY}/ ${WWW}/${FOLDER_DOXY}
	mkdir -p ${WWW}/${FOLDER_DOCB}
	rsync -avz docbook/${FOLDER_DOCB}/ ${WWW}/${FOLDER_DOCB}
	mkdir -p ${WWW}/planning
	rsync -avz planning/ ${WWW}/planning/
