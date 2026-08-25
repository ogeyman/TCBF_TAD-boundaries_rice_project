library(googlesheets4)
gs4_auth()
TAD_groups <- read_sheet("1nYCwVnE4Fm0XRWNQuYxQAS7or4EjZD0o6VKaBO47dUc",sheet="Boundaries_Groups")
TAD_groups <- TAD_groups[,-1]

#meridionalis
a = sum(TAD_groups[,1]=="0") +
sum(TAD_groups[,1]=="1") +
sum(TAD_groups[,1]=="2")+
sum(TAD_groups[,1]=="3")+
sum(TAD_groups[,1]=="4")

#azucena
b = sum(TAD_groups[,2]=="0")+
sum(TAD_groups[,2]=="1")+
sum(TAD_groups[,2]=="2")+
sum(TAD_groups[,2]=="3")+
sum(TAD_groups[,2]=="4")

#rufipogon
c = sum(TAD_groups[,3]=="0")+
sum(TAD_groups[,3]=="1")+
sum(TAD_groups[,3]=="2")+
sum(TAD_groups[,3]=="3")+
sum(TAD_groups[,3]=="4")

#nipponbare
d = sum(TAD_groups[,4]=="0")+
sum(TAD_groups[,4]=="1")+
sum(TAD_groups[,4]=="2")+
sum(TAD_groups[,4]=="3")+
sum(TAD_groups[,4]=="4")

#IR64
e = sum(TAD_groups[,5]=="0")+
sum(TAD_groups[,5]=="1")+
sum(TAD_groups[,5]=="2")+
sum(TAD_groups[,5]=="3")+
sum(TAD_groups[,5]=="4")

a+b+c+d+e
1630*5 # # of rows x # of columns
#make sure line 41 = line 42
