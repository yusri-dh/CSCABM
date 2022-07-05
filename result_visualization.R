# load library ------------------------------------------------------------

setwd("~/Documents/project/agent_4")
library(tidyverse)
library(latex2exp)
library(ggpubr)
# load csv file 1 ---------------------------------------------------------

df = read_csv("./summary_df_mean.csv") %>%
  filter(exp_type != "targeted_therapy_2Gy_24h_30days_mut_05") %>%
  rename(
    CSC_ratio = isStem_function,
    quiescent_ratio = isQuiescent_function,
    n_cells = nrow,
  ) %>%
  mutate(mut_rate = as.double(str_sub(exp_type,-2,-1)) * 0.1) %>%
  mutate(condition=(str_split(exp_type,"_")%>% map_chr(., 1))) 


long_df = gather(df,parameter,value,p_proliferate_mean:n_cells,factor_key=TRUE)

long_df$parameter <- factor(long_df$parameter,labels=c(
  'p_proliferate_mean'=parse(text=TeX("$p_{prol}$")),
  'p_migrate_mean'=parse(text=TeX("$p_{mig}$")),
  'p_symmetry_mean'=parse(text=TeX("$p_{sym}$")),
  'p_dedifferentiation_mean'=parse(text=TeX("$p_{dediff}$")),
  'resistance_mean'=parse(text=TeX('resistance $\\lambda$')),
  'p_capacity_mean'=parse(text=TeX("$p_{cap}$")),
  "CSC_ratio"=parse(text=TeX('$CSC\\,ratio$')),
  "quiescent_ratio"=parse(text=TeX('$quiescent\\,cell\\,ratio$')),
  "n_cells"=parse(text=TeX('$cell\\,number$'))
  ))





# data vis mut 0.05 ------------------------------------------------------------------------------
f <- ggplot(long_df%>% 
              filter(mut_rate==0.5)%>% 
              filter(parameter !="CSC * phantom(.) * number ") %>%
              filter(parameter !="CC * phantom(.) * number "),
            aes(x=step,y=value,color=exp_type))

f1 <- f + geom_smooth(method="gam",alpha=0.5) +
  geom_vline(aes(xintercept = 720),linetype = "longdash")+
  labs(x="time step",y="parameter value")+
  theme_pubr()+
  facet_wrap( ~ parameter,scales = "free_y",
              labeller=label_parsed
              )+
  scale_color_discrete(name="",
                       breaks=c("control_mut_05", "therapy_1Gy_12h_30days_mut_05",
                                "therapy_2Gy_24h_30days_mut_05", "therapy_4Gy_48h_30days_mut_05"),
                       labels=c("Control", "Hyperfractionated radiotherapy", 
                                "Conventional radiotherapy","Hypofractionated radiotherapy")) +
  theme(legend.position = "bottom",plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"))
  
f1
image_path="/home/yusri/Documents/project/latex_agent3/Template for submissions to Scientific Reports/fig"
ggsave("traits_evolution.svg",f1,width = 10,height=6.3,path=image_path)
ggsave("traits_evolution.png",f1,width = 10,height=6.3,dpi=300,path=image_path)

# data vis mutation comparison----------------------------------------------------------------
mut_comparison_df = long_df %>% 
  filter(!exp_type %in% c("therapy_1Gy_12h_30days_mut_05","therapy_4Gy_48h_30days_mut_05"))%>% 
           filter(parameter !="n_CSC") %>%
           filter(parameter !="n_CC") 

f <- ggplot(mut_comparison_df,
                aes(x=step,y=value,color=factor(mut_rate),lty=condition))

f2 <- f + geom_smooth(se=TRUE) +
  geom_vline(aes(xintercept = 720), linetype = "longdash")+
  labs(x="time step",y="parameter value")+
  theme_pubr()+
  facet_wrap( ~parameter,scales = "free_y",
              labeller=label_parsed
  ) +
  scale_color_discrete(name=c("Mutation rate"))+
  scale_linetype_discrete(name="Group", labels=c("Control", "Irradiated"))+
  
  theme(legend.position = "bottom",plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"))
f2
ggsave("mut_comparison.svg",f2,width = 10,height=6.3,dpi=300,path=image_path)
ggsave("mut_comparison.png",f2,width = 10,height=6.3,dpi=300,path=image_path)


# targteted survival ------------------------------------------------------
tgt_df = read_csv("./summary_df_mean.csv") %>%
  filter(exp_type %in% c("targeted_therapy_2Gy_24h_30days_mut_05",
                         "control_mut_05",
                         "therapy_2Gy_24h_30days_mut_05")) %>%
  rename(
    CSC_ratio = isStem_function,
    quiescent_ratio = isQuiescent_function,
    n_cells = nrow,
  ) %>%
  mutate(condition= case_when(
    exp_type=="targeted_therapy_2Gy_24h_30days_mut_05" ~ "targeted radiotherapy",
    exp_type=="control_mut_05" ~ "control",
    exp_type=="therapy_2Gy_24h_30days_mut_05" ~ "conventional radiotherapy"
  )
         )
tgt_df_survival = tgt_df %>% 
  group_by(condition,exp_id) %>%
  summarize(survival_day = max(step))


# targted_vis -------------------------------------------------------------
f <- ggplot(tgt_df,
            aes(x=step,y=n_cells,color=factor(exp_id)))

f3 <- f + geom_line(alpha=0.5) +
  geom_vline(aes(xintercept = 720), linetype = "longdash")+
  scale_y_continuous(trans = "log1p",breaks = c(1,1e1,1e2,1e3,1e4,1e5))+
  facet_wrap( ~condition)+
  theme_pubr()+
  labs(y="log(cell number + 1)")+
  theme(legend.position = "none" ,plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"))
f3

ggsave("targeted_comparison_a.svg",f3,width = 10,height=5,dpi=300,path=image_path)
# targted_vis2 -------------------------------------------------------------
f <- ggplot(tgt_df_survival %>% 
  filter(condition=="targeted radiotherapy"),aes(x=survival_day))
f4 <- f+geom_histogram(binwidth=25) +
  geom_vline(aes(xintercept = 720),linetype = "longdash")+
  scale_x_continuous(name ="Survival day", breaks=seq(400,1200,50))
f4

ggsave("targeted_comparison_b.svg",f4,width = 10,height=6.3,dpi=300,path=image_path)

##
long_df2 = df %>%
  mutate(n_CSC=round(CSC_ratio*n_cells),.before=n_cells) %>%
  mutate(n_CC = n_cells-n_CSC,.before=n_cells) %>% 
  gather(parameter,value,p_proliferate_mean:n_cells,factor_key=TRUE)%>%
  filter(parameter %in% c("n_CSC","n_CC"))%>%
  filter(exp_type %in% c("therapy_2Gy_24h_30days_mut_05",
                         "therapy_1Gy_12h_30days_mut_05",
                         "therapy_4Gy_48h_30days_mut_05"))


f <- ggplot(long_df2  ,aes(x=step,y=value,color = parameter))
f5 =  f + geom_smooth(method="gam",alpha=0.5) +
  geom_vline(aes(xintercept = 720),linetype = "longdash")+
  labs(x="time step",y="cell number")+
  facet_wrap( ~ exp_type)+
  theme_pubr()+
  scale_color_discrete(name="",
                       labels=c("CSC", "CC")) +
  theme(legend.position = "bottom",plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"))

f5

##

grouped_long_df2 = group_by(long_df2%>%filter(step=),exp_type,exp_id,step,parameter,step)
for (i in 1:20) {
  CSC_t0 = long_df2 %>% filter(exp_type == "therapy_2Gy_24h_30days_mut_05" & 
                                 step==0 & 
                                 exp_id==i &
                                 parametyer== "n_CSC")%>% pull(value)
  CSC_t0 = long_df2 %>% filter(exp_type == "therapy_2Gy_24h_30days_mut_05" & 
                                 step==0 & 
                                 exp_id==i &
                                 parametyer== "n_CSC")%>% pull(value)
}