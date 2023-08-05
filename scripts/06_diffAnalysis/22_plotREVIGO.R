######################################################################
### Goal: Edit REVIGO plots 
### Author: Janay Fox
### R script
#####################################################################

## Set up ##
library(ggplot2)
library(scales)
library(dplyr)
library(viridis)
library(ggrepel)
library(patchwork)

#load in data
BN.mfuzz.cluster6.revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
BN.mfuzz.cluster6.revigo.data <- rbind(c("GO:0006575","cellular modified amino acid metabolic process",0.6851865747895725,-1.4626573931674696,-7.375957401785821,5.313865112138396,-3.275537877533078,0.8045735065486667,0.12189312),
                     c("GO:0018126","protein hydroxylation",0.014359122143742097,-4.660596939435513,0.46732149939654,3.635282637998212,-5.976993696766804,0.7625601204623326,0.07550465),
                     c("GO:0018208","peptidyl-proline modification",0.20196901915454343,-4.437312255783912,-2.4199390646681382,4.783346067476675,-5.06095941946776,0.6636365526335498,0.48991279),
                     c("GO:0019471","4-hydroxyproline metabolic process",0.00013969959000166042,5.223544575637416,-2.965993957712103,1.6334684555795864,-6.709255333356367,0.8080486425796953,0),
                     c("GO:0019511","peptidyl-proline hydroxylation",0.008737876736532427,-5.962441729143965,-1.722079985237315,3.419625360887743,-6.5960949160448035,0.623059453497824,0.24916865),
                     c("GO:0071456","cellular response to hypoxia",0.013477684254445906,1.9735185641624111,5.149955091525602,3.607776603741693,-2.1106418862833567,0.6559244350235264,0.00436528),
                     c("GO:1901605","alpha-amino acid metabolic process",3.5585210895822956,1.495944422195505,-3.9409029151515487,6.029324108830442,-2.904129865219576,0.66522545048598,0.32679993));

BN.mfuzz.cluster5.revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
BN.mfuzz.cluster5.revigo.data <- rbind(c("GO:0002252","immune effector process",0.060752690747150666,-5.998094700820335,-4.598450761702338,4.26164345331189,-300,0.6216655149410969,0.68371319),
                     c("GO:0002376","immune system process",0.5358776487042264,-0.2775107065212202,6.449903227701661,5.207122497650964,-4.375068085808162,1,-0),
                     c("GO:0002682","regulation of immune system process",0.29494907484279137,0.38160185579835487,-6.426453588376232,4.947806094714132,-2.3318140635124465,0.8503288591443282,0.12346543),
                     c("GO:0006952","defense response",1.0410213876080876,-7.070536298917201,1.6169622180610175,5.495515198600499,-3.0263368029776636,0.6752063450058148,0.26469029),
                     c("GO:0006956","complement activation",0.02405161274528587,-5.15145922938502,-2.4184092854987003,3.859258417467307,-300,0.24511752745717946,0.68120383),
                     c("GO:0006959","humoral immune response",0.040652580690483185,-6.3910942510337945,-1.7782713508850148,4.087177811761157,-300,0.5212241731017286,-0),
                     c("GO:0009617","response to bacterium",0.14606922606959327,-5.82063423423205,2.629817287853001,4.642622776409274,-2.8459599313518695,0.5644667023925283,0.29520673),
                     c("GO:0072376","protein activation cascade",0.0026675969328888494,3.8045675074541405,2.436255274690121,2.904715545278681,-300,1,0),
                     c("GO:2000257","regulation of protein activation cascade",0.00020954938500249067,3.6730697752147288,-3.2658876857394574,1.806179973983887,-2.5539372080910536,0.896905571765112,-0));


BA.mfuzz.cluster5.revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
BA.mfuzz.cluster5.revigo.data <- rbind(c("GO:0001666","response to hypoxia",0.03421642100826383,-4.0564464713848185,5.205032264723029,4.012330955580147,-1.9071272350098294,0.7168140071821428,0.13752985),
                     c("GO:0002412","antigen transcytosis by M cells in mucosal-associated lymphoid tissue",1.3304722857300992E-05,-6.0776870834810905,1.269066500353262,0.6989700043360189,-2.4567900040092248,0.8291741136955408,-0),
                     c("GO:0018126","protein hydroxylation",0.014359122143742097,0.1285268881145624,-6.58437813043557,3.635282637998212,-1.698099181216368,0.9320867795910414,0.24916865),
                     c("GO:0018208","peptidyl-proline modification",0.20196901915454343,1.0901376541411922,-6.284790046967555,4.783346067476675,-1.3344787881880034,0.8955853097897764,0.48991279),
                     c("GO:0019471","4-hydroxyproline metabolic process",0.00013969959000166042,-3.9188251835338996,-7.11707987870439,1.6334684555795864,-1.9272992614619782,0.9749626429166117,-0),
                     c("GO:0019511","peptidyl-proline hydroxylation",0.008737876736532427,1.6548523068076866,-6.256654004212935,3.419625360887743,-1.9071272350098294,0.867255815242152,0.07386162),
                     c("GO:0032364","intracellular oxygen homeostasis",0.0007716739257234576,0.7042894830821558,-0.7491130199456139,2.367355921026019,-2.719521279912308,1,-0),
                     c("GO:0045056","transcytosis",0.0023283265000276737,-7.038519184452494,-0.0193563275157307,2.8457180179666586,-1.834765941803245,0.8673164785002478,0.35633827),
                     c("GO:0051342","regulation of cyclic-nucleotide phosphodiesterase activity",0.0004756438421485105,5.56266520482584,-2.342156089961538,2.1583624920952498,-2.321121859120142,0.938259661606441,0.42930321),
                     c("GO:0051344","negative regulation of cyclic-nucleotide phosphodiesterase activity",0.0001596566742876119,5.404183837073713,-2.941804983460242,1.6901960800285136,-2.719521279912308,0.9393948247648263,0),
                     c("GO:0055008","cardiac muscle tissue morphogenesis",0.008621460411531045,-6.2569946270004895,-1.8258014961697053,3.4138025167693513,-1.906509695903187,0.6803546744097171,0.63730513),
                     c("GO:0060343","trabecula formation",0.002624356583602621,-6.146905667207981,-3.5975843796489153,2.8976270912904414,-1.906509695903187,0.7532080836154397,0.49804856),
                     c("GO:0060711","labyrinthine layer development",0.00529860587792012,-6.611596192845193,-1.5717684519346242,3.2024883170600935,-2.2162621937396665,0.7604810534603007,0.28466572),
                     c("GO:0070482","response to oxygen levels",0.038207837865454126,-4.4397523317077185,5.289418276049121,4.060244426898225,-1.7862744074685424,0.7848039330322394,0.65312022),
                     c("GO:0071731","response to nitric oxide",0.0014236053457312061,-2.2464802711517016,6.132690705384817,2.6324572921847245,-1.906509695903187,0.900363156703191,0.17110639),
                     c("GO:0099159","regulation of modification of postsynaptic structure",0.0007151288535799284,2.24252873665915,6.508583064952034,2.3344537511509307,-2.628745250794901,0.9625634410252221,0.08962698),
                     c("GO:0140252","regulation protein catabolic process at postsynapse",0.0003791846014330783,1.3076389718485109,3.6345422073414815,2.060697840353612,-2.586088988152483,0.9644350427350025,0.09591616),
                     c("GO:1905289","regulation of CAMKK-AMPK signaling cascade",0.00038251078214740356,5.181074105828184,2.783842343395377,2.0644579892269186,-1.7981922431626012,0.9206081785744782,0.10402292),
                     c("GO:1905290","negative regulation of CAMKK-AMPK signaling cascade",6.652361428650496E-06,5.369031897330604,2.3319599157407804,0.47712125471966244,-1.906509695903187,0.924652616247932,0.63661272));


## Parse data ##
revigo.dat.parse <- function(revigo.data, revigo.names){
  one.data <- data.frame(revigo.data)
  names(one.data) <- revigo.names
  one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ]
  one.data$plot_X <- as.numeric( as.character(one.data$plot_X) )
  one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) )
  one.data$log_size <- as.numeric( as.character(one.data$log_size) )
  one.data$value <- as.numeric( as.character(one.data$value) )
  one.data$frequency <- as.numeric( as.character(one.data$frequency) )
  one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) )
  one.data$dispensability <- as.numeric( as.character(one.data$dispensability) )
  return(one.data)
}

BN.mfuzz.cluster6.one.data <- revigo.dat.parse(BN.mfuzz.cluster6.revigo.data, BN.mfuzz.cluster6.revigo.names)
BN.mfuzz.cluster5.one.data <- revigo.dat.parse(BN.mfuzz.cluster5.revigo.data, BN.mfuzz.cluster5.revigo.names)

BA.mfuzz.cluster5.one.data <- revigo.dat.parse(BA.mfuzz.cluster5.revigo.data, BA.mfuzz.cluster5.revigo.names)

#convert p value and size columns 
BN.mfuzz.cluster6.one.data$value  <- as.numeric(BN.mfuzz.cluster6.one.data$value %>% lapply(function(x) 10^x))
BN.mfuzz.cluster5.one.data$value  <- as.numeric(BN.mfuzz.cluster5.one.data$value %>% lapply(function(x) 10^x))

BA.mfuzz.cluster5.one.data$value  <- as.numeric(BA.mfuzz.cluster5.one.data$value %>% lapply(function(x) 10^x))

## Plot ##  
plot.revigo <- function(one.data){
  ex <- one.data [ one.data$dispensability < 0.35, ]
  one.x_range = max(one.data$plot_X) - min(one.data$plot_X)
  one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y)
  p1 <- ggplot( data = one.data ) + geom_point( aes( plot_X, plot_Y, color = value, size = log_size), alpha = I(0.8) ) + 
    scale_colour_viridis(option = "viridis", limits = c(0, 0.05), name = "q-value" , direction = -1) + 
    geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + 
    theme_bw() + geom_text_repel(data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 4, 
                                  min.segment.length = Inf, max.overlaps = 20) +
    labs (y = "Semantic space y", x = "Semantic space x") + theme(legend.key = element_blank()) + 
    xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10) + 
    ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10) + 
    scale_size_continuous(breaks = c(1, 3, 5, 7), name = "Log Size", limits = c(0, 7)) + 
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), legend.text = element_text(size=11), 
          legend.position = "right")
  return(p1)
}

BN.mfuzz.cluster6.plot <- plot.revigo(BN.mfuzz.cluster6.one.data)
BN.mfuzz.cluster5.plot <- plot.revigo(BN.mfuzz.cluster5.one.data)

BA.mfuzz.cluster5.plot <- plot.revigo(BA.mfuzz.cluster5.one.data)

#save plots 
ggsave("BN_StFvsSwF_up_REVIGO.tiff", plot = BN.StFvsSwF.up.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BN_StFvsSwF_down_REVIGO.tiff", plot = BN.StFvsSwF.down.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BN_SwFvsSwP_down_REVIGO.tiff", plot = BN.SwFvsSwP.down.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BN_vennGenes_REVIGO.tiff", plot = BN.venngenes.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BN_mfuzz_cluster6_0.35_REVIGO.tiff", plot = BN.mfuzz.cluster6.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BN_mfuzz_cluster5_0.35_REVIGO.tiff", plot = BN.mfuzz.cluster5.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BN_30perc_REVIGO.tiff", plot = BN.30perc.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BN_sameDir_up_0.35_REVIGO.tiff", plot = BN.sameDir.up.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)

ggsave("BA_StFvsSwF_up_REVIGO.tiff", plot = BA.StFvsSwF.up.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BA_StFvsSwF_down_REVIGO.tiff", plot = BA.StFvsSwF.down.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BA_SwFvsSwP_up_REVIGO.tiff", plot = BA.SwFvsSwP.up.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BA_SwFvsSwP_down_REVIGO.tiff", plot = BA.SwFvsSwP.down.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BA_mfuzz_cluster5_0.35_REVIGO.tiff", plot = BA.mfuzz.cluster5.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
ggsave("BA_30perc_REVIGO.tiff", plot = BA.30perc.plot, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)

# Panel BN plots 
panel.plot <- (BN.mfuzz.cluster5.plot +  BN.mfuzz.cluster6.plot) / BA.mfuzz.cluster5.plot +
                  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
panel.plot

ggsave("REVIGO_panel.tiff", plot = panel.plot, device = "tiff", width = 11, height = 7, units = "in", dpi = 600)
