# Set the working in directory to control where plots are saved to
setwd("~/Desktop/")
require(ggplot2)
# use to make multipanel ggplots
require(patchwork)
require(deSolve)

beta <- function(t){
  beta = 2
  if(t>10){
    beta=0.01
  }
  return(beta)
}

sigma <- 1/5
gamma <- 1/5
h <- 1/100
rho <- 0.3
eta <- 1/2
delta <- 1/5
T<-60

epi<-function(t, y, parameters){
  E = y[1]
  A = y[2]
  I = y[3]
  H = y[4]
  C = y[5]
  
dE= beta(t)*(A + I) - sigma*E
dA= rho*sigma*E - gamma*A
dI = (1-rho)*sigma*E - (1-h)*gamma*I - h*eta*I
dH = h*eta*I - delta*H
dC = (1-rho)*sigma*E
return(list(c(dE,dA,dI,dH,dC)))
}

output <- ode(y = c(E=10,A=5,I=5,H=0, C=5), parms = NULL, times = seq(0,T), func=epi)
output <- data.frame(output)
incidence <- diff(c(0,output$C))
output <- data.frame(output, incidence = incidence)

beta <- function(t){
  beta = 2
  if(t>10){
    beta=0.12
  }
  return(beta)
}

output2 <- ode(y = c(E=10,A=5,I=5,H=0, C=5), parms = NULL, times = seq(0,T), func=epi)
output2 <- data.frame(output2)
incidence <- diff(c(0,output2$C))
output2 <- data.frame(output2, incidence = incidence)

times = seq(0,T, .1)
i = min(which(times>10))
strict = rep(-1, length(times))
mild = strict
strict[i:length(times)]=0
mild[1:(i-1)] = 0
measures = data.frame(times, strict = strict, mild = mild)

gha=ggplot()+
  geom_ribbon(data=output, aes(x=time, ymax=I*max(H)/max(I), ymin=0), fill = "grey50", alpha=0.3)+
  geom_line(data=output, aes(x=time, y=H), color = "purple")+
  geom_ribbon(data = measures, aes(x=times,ymin = -.5, ymax = .5*strict), fill = "red")+
  geom_ribbon(data = measures, aes(x=times,ymin = -.5, ymax = .5*mild), fill = "green")+
  ggtitle("Very effective strict NPIs")+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("hospital occupany")+
  scale_x_continuous(breaks=seq(min(output$time),max(output$time),3))+
  xlab("days")

ghb = ggplot()+
  geom_ribbon(data=output2, aes(x=time, ymax=I*max(H)/max(I), ymin=0), fill = "grey50", alpha=0.3)+
  geom_line(data=output2, aes(x=time, y=H), color = "purple")+
  geom_ribbon(data = measures, aes(x=times,ymin = -.58, ymax = .58*strict), fill = "red")+
  geom_ribbon(data = measures, aes(x=times,ymin = -.58, ymax = .58*mild), fill = "green")+
  ggtitle("Less effective strict NPIs")+
  scale_x_continuous(breaks=seq(min(output$time),max(output$time),3))+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("hospital occupancy")+
  xlab("days")


gca=ggplot()+
  geom_ribbon(data=output, aes(x=time, ymax=I*max(incidence)/max(I), ymin=0), fill = "grey50", alpha=0.3)+
  geom_line(data=output, aes(x=time, y=incidence), color = "dodgerblue")+
  geom_ribbon(data = measures, aes(x=times,ymin = -20, ymax = 20*strict), fill = "red")+
  geom_ribbon(data = measures, aes(x=times,ymin = -20, ymax = 20*mild), fill = "green")+
  ylab("new cases")+
  scale_x_continuous(breaks=seq(min(output$time),max(output$time),3))+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Very effective strict NPIs")+
  xlab("days")

gcb=ggplot()+
  geom_ribbon(data=output2, aes(x=time, ymax=I*max(incidence)/max(I), ymin=0), fill = "grey50", alpha=0.3)+  
  geom_line(data=output2, aes(x=time, y=incidence), color = "dodgerblue")+
  geom_ribbon(data = measures, aes(x=times,ymin = -10, ymax = 10*strict), fill = "red")+
  geom_ribbon(data = measures, aes(x=times,ymin = -10, ymax = 10*mild), fill = "green")+
  ylab("new cases")+
  scale_x_continuous(breaks=seq(min(output$time),max(output$time),3))+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Less effective strict NPIs")+
  xlab("days")



###

require(dplyr)
require(chron)
require(ggplot2)
require(patchwork)
data = read.csv("https://raw.githubusercontent.com/ahurford/NL-public-COVID-data/main/NL_prov_stats.csv")
dates <- data$date_of_update
data$date_of_update <- format(as.POSIXct(dates,format='%Y/%m/%d %H:%M:%S'),format='%Y/%m/%d')
data = filter(data,date_of_update >= "2021/02/08")%>%filter(date_of_update < "2021/04/01")%>%select(date_of_update, currently_hospitalized, active_cases, new_provincial_cases)
dates = seq(as.Date("2021/02/08"), as.Date("2021/03/31"), "days")
data = data.frame(data, dates)
# on Avalon
dates2 = c(dates[1:5], dates[5]+0.1, dates[6:34], dates[34]+0.1, dates[35:48], dates[48]+.1, dates[49:length(dates)])
level_5 = rep(-1,length(dates2))
level_5[6:35] = 0
level_2 = rep(-1,length(dates2))
level_2[c(1:5,51:length(dates2))] = 0
level_4 = rep(-1,length(dates2))
level_4[36:50] = 0
alert_levels = data.frame(dates2, level_5=level_5,level_2=level_2,level_4=level_4)
# March 13 - level 4 (on Avalon)
# Feb 12 - level 5 (whole province - although some restrictions prior)
# March 27 - Alert level 2.

g1=ggplot(data=data, aes(x = dates))+
  geom_ribbon(aes(ymax=11*active_cases/max(active_cases),ymin=0, group=1), fill = "grey70", alpha =0.5)+
  geom_line(aes(y=currently_hospitalized, group=1), color = "purple")+
  xlab("")+
  geom_ribbon(data = alert_levels, aes(x=dates2,ymin = -.5, ymax = .5*level_2), fill = "green")+
  geom_ribbon(data = alert_levels, aes(x=dates2,ymin = -.5, ymax = .5*level_5), fill = "red")+
  geom_ribbon(data = alert_levels, aes(x=dates2,ymin = -.5, ymax = .5*level_4), fill = "orange")+
  scale_x_date(date_labels="%d %b",date_breaks  ="3 day")+
  ggtitle("Alpha variant outbreak in Nfld.")+
  xlab("2021")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")+
  ylab("currently hospitalized")



g2 = ggplot(data=data, aes(x = dates))+
  geom_ribbon(aes(ymax=max(new_provincial_cases)*active_cases/max(active_cases),ymin=0, group=1), fill = "grey70", alpha =0.5)+
  geom_line(aes(y=new_provincial_cases, group=1), color = "dodgerblue")+
  xlab("")+
  geom_ribbon(data = alert_levels, aes(x=dates2,ymin = -5, ymax = 5*level_2), fill = "green")+
  geom_ribbon(data = alert_levels, aes(x=dates2,ymin = -5, ymax = 5*level_5), fill = "red")+
  geom_ribbon(data = alert_levels, aes(x=dates2,ymin = -5, ymax = 5*level_4), fill = "orange")+
  ylab("new cases")+
  ggtitle("Alpha variant outbreak in Nfld.")+
  xlab("2021")+
  scale_x_date(date_labels="%d %b",date_breaks  ="3 day")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")

g3=(g2+g1)
g1=(gca+gha)
g2=(gcb+ghb)
g1/g2/g3+plot_annotation(tag_levels = 'A')
ggsave("figures/scenarios.png", units = "cm", height = 15, width=18)

