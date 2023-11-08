###########################    FIGURE 2    #####################################
################################################################################
################################################################################
######  This script is used to extract slope-coefficients to reproduce the boxplots
###### showing the bias of the different scenarios.
###### Simulation results were executed on a compute server (cf. readme).  
#   - Binary and Continuous data
#   - Sample sizes of N = 2000, 6000, and 20K
#   - Different values of alpha: 2, 3, and 4
#   - Four distinct methods: Naive, Log/Residual, Linear/Substitution, and Log/Substitution
#
###### The extracted coefficients are used in constructing graphs, which graphically 
###### represent the bias of these coefficients over time, under different scenarios.
###### Relative bias boxplot (corresponding to Figure 2) are plotted in file
###### "Figure2.pdf"
################################################################################
################################################################################

##### Loading packages #####
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(here)

##### Define directory #####
directory <- here()
setwd(directory)

##### Loading simulated data  #####
#For logistic/substitution
getwd()
setwd("./log_sub")
pattern <- c("simu_boots_binaire_2000_")
files <- list.files(pattern=pattern)
pattern <- "simu_boots_binaire_6000_"
files1 <- list.files(pattern=pattern)
pattern <- "simu_boots_binaire_20000_"
files2 <- list.files(pattern=pattern)
# For linear substitution
setwd("../linear_sub")
pattern <- "simu_boots_2000_LINEAR"
files3 <- list.files(pattern=pattern)
pattern <- "simu_boots_6000_LINEAR"
files4 <- list.files(pattern=pattern)
pattern <- "simu_boots_20000_LINEAR"
files5 <- list.files(pattern=pattern)
# For logistic/Residual
setwd("../log_res")
pattern <- "simu_boots_2sri_2000_"
files6 <- list.files(pattern=pattern)
pattern <- "simu_boots_2sri_6000_"
files7 <- list.files(pattern=pattern)
pattern <- "simu_boots_2sri_20000_"
files8 <- list.files(pattern=pattern)
files <- c(files,files1,files2,files3,files4,files5,files6,files7,files8)
# For CONT/LOG-SUB
setwd("../cont_log_sub")
pattern <- "simu_boots_2000.R"
files_c_1 <- list.files(pattern=pattern)
pattern <- "simu_boots_6000_05.R"
files_c_2 <- list.files(pattern=pattern)
pattern <- "simu_boots_05_20000_"
files_c_3 <- list.files(pattern = pattern)
pattern <- "simu_boots_2000_4"
files_c_4 <- list.files(pattern=pattern)
pattern <- "simu_boots_6000.R"
files_c_5 <- list.files(pattern=pattern)
pattern <- "simu_boots_20000_"
files_c_6 <- list.files(pattern = pattern)
files_c <- c(files_c_1,files_c_2,files_c_3,files_c_4,files_c_5,files_c_6)

##### Combine results : #####
coefN_inter <- NULL
coefIV_inter <- NULL
setwd("../log_sub")
r <- 0
for(f in files[1:4500]){
  r <- r+1
  load(f)
  coefN_inter <- append(coefN_inter,as.numeric(res[seq(5,length(res),25)]))
  coefIV_inter <- append(coefIV_inter,as.numeric(res[seq(10,length(res),25)]))
}
setwd("../linear_sub")
for(f in files[4501:9000]){
  r <- r+1
  load(f)
  coefIV_inter <- append(coefIV_inter,as.numeric(res[seq(10,length(res),26)]))
}
setwd("../log_res")
for(f in files[9001:13500]){
  r <- r+1
  load(f)
  coefIV_inter <- append(coefIV_inter,as.numeric(res[seq(6,length(res),21)]))
}

r <- 0
coefIV_inter_c <- NULL
coefN_inter_c <- NULL
setwd("../cont_log_sub")
for(f in files_c){
  r <- r+1
  load(f)
coefN_inter_c <- append(coefN_inter_c,as.numeric(res[seq(5,length(res),26)]))
coefIV_inter_c <- append(coefIV_inter_c,as.numeric(res[seq(10,length(res),26)]))
}

##### Graphics generation : ####
alpha <- c(rep("Naive",1500), rep("Logistic/substitution", 1500), rep("Linear/substitution", 1500), rep("Logistic/residual", 1500))
alpha <- factor(alpha , levels=c("Naive","Logistic/residual", "Logistic/substitution","Linear/substitution"))
names <- c(rep("2000", 500) , rep("6000", 500) , rep("20K", 500))
names <- factor(names , levels=c("2000", "6000","20K"))

value <- c((coefN_inter[1:500]/rep(1,500)-1),
           (coefN_inter[1501:2000]/rep(1,500)-1),
           (coefN_inter[3001:3500]/rep(1,500)-1),
           (coefIV_inter[1:500]/rep(1,500)-1),
           (coefIV_inter[1501:2000]/rep(1,500)-1),
           (coefIV_inter[3001:3500]/rep(1,500)-1),
           (coefIV_inter[4501:5000]/rep(1,500)-1),
           (coefIV_inter[6001:6500]/rep(1,500)-1),
           (coefIV_inter[7501:8000]/rep(1,500)-1),
           (coefIV_inter[9001:9500]/rep(1,500)-1),
           (coefIV_inter[10501:11000]/rep(1,500)-1),
           (coefIV_inter[12001:12500]/rep(1,500)-1))

data <- data.frame(names,value,alpha)


p2 <- ggplot(data, aes(x=alpha,y=value, fill=names))+
  geom_boxplot()+
  ylab("Relative_bias")+
  ylim(-0.5,1.5)+
  scale_fill_manual(values=c("darkslategrey","darkslategray4","grey"))+
  labs(fill = "",title= expression("Binary ; " ~alpha~"= 2 ; R² = 14.3"),x ="" , y = "Relative bias")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),
        legend.title = element_text(color = "black", size = 25),
        legend.text = element_text(color = "black", size = 20),
        legend.key.size = unit(0.9, "cm"),
        legend.key.width = unit(0.9,"cm"), 
        legend.position ="none",
        plot.title = element_text(color="black", size=20, face="bold", ),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15))


value <- c((coefN_inter[501:1000]/rep(1,500)-1),
           (coefN_inter[2001:2500]/rep(1,500)-1),
           (coefN_inter[3501:4000]/rep(1,500)-1),
           (coefIV_inter[501:1000]/rep(1,500)-1),
           (coefIV_inter[2001:2500]/rep(1,500)-1),
           (coefIV_inter[3501:4000]/rep(1,500)-1),
           (coefIV_inter[5001:5500]/rep(1,500)-1),
           (coefIV_inter[6501:7000]/rep(1,500)-1),
           (coefIV_inter[8001:8500]/rep(1,500)-1),
           (coefIV_inter[9501:10000]/rep(1,500)-1),
           (coefIV_inter[11001:11500]/rep(1,500)-1),
           (coefIV_inter[12501:13000]/rep(1,500)-1))

data <- data.frame(names,value,alpha)

p3 <- ggplot(data, aes(x=alpha,y=value, fill=names))+
  geom_boxplot()+
  ylab("Relative_bias")+
  ylim(-0.5,1.5)+
  scale_fill_manual(values=c("darkslategrey","darkslategray4","grey"))+
  labs(fill = "",title= expression("Binary ; " ~alpha~"= 3 ; R² = 35.2"),x ="" , y = "Relative bias")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),
        legend.title = element_text(color = "black", size = 25),
        legend.text = element_text(color = "black", size = 20),
        legend.key.size = unit(0.9, "cm"),
        legend.key.width = unit(0.9,"cm"), 
        legend.position ="none",
        plot.title = element_text(color="black", size=20, face="bold", ),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15))

value <- c((coefN_inter[1001:1500]/rep(1,500)-1),
           (coefN_inter[2501:3000]/rep(1,500)-1),
           (coefN_inter[4001:4500]/rep(1,500)-1),
           (coefIV_inter[1001:1500]/rep(1,500)-1),
           (coefIV_inter[2501:3000]/rep(1,500)-1),
           (coefIV_inter[4001:4500]/rep(1,500)-1),
           (coefIV_inter[5501:6000]/rep(1,500)-1),
           (coefIV_inter[7001:7500]/rep(1,500)-1),
           (coefIV_inter[8501:9000]/rep(1,500)-1),
           (coefIV_inter[10001:10500]/rep(1,500)-1),
           (coefIV_inter[11501:12000]/rep(1,500)-1),
           (coefIV_inter[13001:13500]/rep(1,500)-1))

data <- data.frame(names,value,alpha)

p4 <- ggplot(data, aes(x=alpha,y=value, fill=names))+
  geom_boxplot()+
  ylab("Relative_bias")+
  ylim(-0.5,1.5)+
  scale_fill_manual(values=c("darkslategrey","darkslategray4","grey"))+
  labs(fill = "",title= expression("Binary ; " ~alpha~"= 4 ; R² = 58.6"),x ="" , y = "Relative bias")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),
        legend.title = element_text(color = "black", size = 25),
        legend.text = element_text(color = "black", size = 20),
        legend.key.size = unit(0.9, "cm"),
        legend.key.width = unit(0.9,"cm"), 
        legend.position ="none",
        plot.title = element_text(color="black", size=20, face="bold", ),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15))

leg <- ggplot(data, aes(x=alpha,y=value, fill=names))+
  geom_boxplot()+
  ylab("Relative_bias")+
  ylim(-0.5,1.5)+
  scale_fill_manual(values=c("darkslategrey","darkslategray4","grey"))+
  labs(fill = "",title= expression("Binary ; " ~alpha~"= 4 ; R² = 58.6"),x ="" , y = "Relative bias")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),
        legend.title = element_text(color = "black", size = 25),
        legend.text = element_text(color = "black", size = 20),
        legend.key.size = unit(0.9, "cm"),
        legend.key.width = unit(0.9,"cm"), 
        legend.position ="bottom",
        plot.title = element_text(color="black", size=20, face="bold", ),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15))

L <- get_legend(leg)

# Graphic for continuous : 

alpha <- c(rep("Naive",1500), rep("Logistic/Substitution", 1500))
alpha <- factor(alpha , levels=c("Naive", "Logistic/Substitution"))
names <- c(rep("2000", 500) , rep("6000", 500) , rep("20K", 500))
names <- factor(names , levels=c("2000", "6000","20K"))
value <- c((coefN_inter_c[1:500]/rep(1,500)-1),
           (coefN_inter_c[501:600]/rep(1,500)-1),
           (coefN_inter_c[601:1100]/rep(1,500)-1),
           (coefIV_inter_c[1:500]/rep(1,500)-1),
           (coefIV_inter_c[501:600]/rep(1,500)-1),
           (coefIV_inter_c[601:1100]/rep(1,500)-1))
data <- data.frame(names,value,alpha)


p5 <- ggplot(data, aes(x=alpha,y=value, fill=names))+
  geom_boxplot()+
  ylim(-0.3,0.5)+
  scale_fill_manual(values=c("darkslategrey","darkslategray4","grey"))+
  labs(fill = "",title= expression("Continuous ; " ~alpha~"= 0.5"),x ="" , y = "Relative bias")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),
        legend.title = element_text(color = "black", size = 25),
        legend.text = element_text(color = "black", size = 20),
        legend.key.size = unit(0.9, "cm"),
        legend.key.width = unit(0.9,"cm"),
        legend.position = "none",
        plot.title = element_text(color="black", size=20, face="bold", ),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15))


value <- c((coefN_inter_c[1101:1600]/rep(1,500)-1),
           (coefN_inter_c[1601:2100]/rep(1,500)-1),
           (coefN_inter_c[2101:2600]/rep(1,500)-1),
           (coefIV_inter_c[1101:1600]/rep(1,500)-1),
           (coefIV_inter_c[1601:2100]/rep(1,500)-1),
           (coefIV_inter_c[2101:2600]/rep(1,500)-1))
data <- data.frame(names,value,alpha)


p6 <- ggplot(data, aes(x=alpha,y=value, fill=names))+
  geom_boxplot()+
  ylab("Relative_bias")+
  ylim(-0.3,0.5)+
  scale_fill_manual(values=c("darkslategrey","darkslategray4","grey"))+
  labs(fill = "",title= expression("Continuous ; " ~alpha~"= 1"),x ="" , y = "Relative bias")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),
        legend.title = element_text(color = "black", size = 25),
        legend.text = element_text(color = "black", size = 20),
        legend.key.size = unit(0.9, "cm"),
        legend.key.width = unit(0.9,"cm"), 
        legend.position = "none",
        plot.title = element_text(color="black", size=20, face="bold", ),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15))

##### Final graph : #####
setwd(directory)
setwd("Figure")
pdf(file ="Figure2.pdf",width=14,height=6)
grid.arrange(p5,p2,p6,p3,L,p4,ncol=2, nrow=3)
dev.off()


