# update the following line according to your installation
SPSK=/path/to/superseeker

all: Sc1toSc9.vcf G1_subclone.vcf G3_subclone.vcf G4_subclone.vcf

Sc1toSc9.vcf: C1toC9.input.vcf
	$(SPSK) $< > $@

G1_subclone.vcf: variantClustersInG1.input.vcf
	$(SPSK) $< > $@

G3_subclone.vcf: variantClustersInG3.input.vcf
	$(SPSK) $< > $@
	
G4_subclone.vcf: variantClustersInG4.input.vcf
	$(SPSK) $< > $@

