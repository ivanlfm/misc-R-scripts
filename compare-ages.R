#fetch ages of equivalent nodes in different trees

library(phytools)
setwd("C:/a/")

#reads trees
iqtree <- read.nexus("iqtree.nex")
mb <-read.nexus("mrbayes.nex")
beast<-read.nexus("beast-exp.nex")

#defines pairs of taxa to get ages
tips1 <- c("M_swainsoni_090",
"M_swainsoni_090",
"M_swainsoni_090",
"M_beta_227",
"M_acuta_004_058",
"M_furcata_037",
"M_furcata_037",
"M_cyanospina_138",
"M_sagittata_187",
"M_banksi_784750",
"M_forcipata_242_DR_785054",
"M_horrida_042",
"M_horrida_042",
"M_horrida_042",
"M_horrida_042",
"M_clypeata_139_169",
"M_vigorsi_140",
"M_vigorsi_140",
"M_pungens_069_070",
"M_aureola_171",
"M_cornuta_145_199",
"M_clypeata_139_169",
"M_crassispina_014",
"M_crassispina_014",
"M_crassispina_014",
"M_crassispina_014",
"M_furva_208",
"M_nigrichelis_120",
"M_nigrichelis_120",
"M_miles_142",
"M_macfarlanei_054_055",
"M_miles_142",
"M_furcula",
"M_furcula",
"M_similis_243_755496",
"M_similis_243_755496",
"M_lucasi_126",
"M_lucasi_126",
"M_shealsi_160",
"M_cf_crassa_213",
"M_cf_crassa_213",
"M_PHMsp01_4265",
"M_PHMsp01_4265",
"M_decorata",
"M_decorata",
"M_duodecimspinosa_218",
"Pronous_peje_JCG")

tips2 <- c("M_duodecimspinosa_218",
"M_triangularispinosa_B_093",
"M_sexspinosa_084_190",
"M_perfida_026",
"M_flaveola_011",
"M_coca_212",
"M_reimoseri_072",
"M_anchicaya_005",
"M_banksi_784750",
"M_militaris_244_DR_784363",
"M_spinulata_205",
"M_spinulata_205",
"M_gracilis_FAPDNA015_00000988A",
"M_horrida_122",
"M_duodecimspinosa_218",
"M_schreibersi_148",
"M_schreibersi_148",
"M_pungens_069_070",
"M_cucharas",
"M_picta_168",
"M_woytkowskii_4262",
"M_duodecimspinosa_218",
"M_duodecimspinosa_218",
"M_sanctispiritus_041",
"M_macfarlanei_054_055",
"M_digitata_016_017",
"M_digitata_016_017",
"M_digitata_016_017",
"M_reali_207",
"M_kirbyi_167",
"M_kirbyi_167",
"M_fissispina_033",
"M_fissispina_033",
"M_cubana_784820",
"M_cubana_784820",
"M_bimucronata_123",
"M_shealsi_160",
"M_duodecimspinosa_218",
"M_pupa_4264",
"M_pupa_4264",
"M_fidelis_192",
"M_gaujoni_039",
"M_pilaton_IFM1712",
"M_lepidoptera_130",
"M_duodecimspinosa_218",
"M_plana_062",
"M_duodecimspinosa_218")

#finds tree heights
THiqtree <- node.depth.edgelength(iqtree)[1]
THmb <- node.depth.edgelength(mb)[1]
THbeast <- node.depth.edgelength(beast)[1]

ages <- matrix(NA, nrow = length(tips1), ncol = 5)
colnames(ages)<-c("terminal1", "terminal2","iqtree", "mrbayes", "beast")
for(i in 1:length(tips[,1])) {
	ages[i,1] <- tips1[i]
	ages[i,2] <- tips2[i]
	ages[i,3] <- THiqtree-findMRCA(iqtree,tips=c(tips1[i], tips2[i]), type="height")
	ages[i,4] <- THmb-findMRCA(mb,tips=c(tips1[i], tips2[i]), type="height")
	ages[i,5] <- THbeast-findMRCA(beast,tips=c(tips1[i], tips2[i]), type="height")
	}