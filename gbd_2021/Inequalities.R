# 健康不公平的可视化
# 2022-12-08

# 加载需要的包 ------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(car)# 异方差诊断
library(MASS)# 稳健回归
# install.packages("devtools")
# devtools::install_github("nicolash2/ggbrace")
library(ggbrace)# 画括号
library(mgcv)# 提供洛伦兹曲线拟合 (样条函数等)
library(splines)# 拟合样条函数
library(broom)

# 数据准备--------------------------------------------------------------------
# 读取并合并三个数据
setwd("C:\\Users\\35111\\Desktop\\健康不公平的可视化")
# 疾病负担数据
burden <- read.csv("数据\\all causes.csv",header = T)
burden <- burden %>%
  filter(age!="Age-standardized")
# sdi 数据
sdi <- read.csv("数据\\各国各年份SDI(已校对).csv",header = T,check.names = F)
sdi <- sdi %>% # 宽数据转为长数据
  pivot_longer(cols = `1990`:`2019`,names_to = "year") %>%
  rename(sdi=value) %>%
  dplyr::select(location,year,sdi)
sdi$year <- as.integer(sdi$year)
# 匹配 sdi 与疾病负担, 生成 data
data <- left_join(burden,sdi,by=c("location","year"))
# 注意一个国家 Georgia
# 读取人口数据
dir <- "C:\\Users\\35111\\Desktop\\健康不公平的可视化\\数据\\GBD人口数据"
files <- list.files(dir,
                    pattern = ".CSV",
                    full.names = T)
pop <- map_dfr(files,read.csv)
pop <- pop %>%
  dplyr::select(location_name,sex_name,age_group_name,year_id,metric_name,val)
colnames(pop) <- c("location","sex","age","year","metric","val")
pop$sex[pop$sex=="male"] <- "Male"
pop$sex[pop$sex=="female"] <- "Female"
pop$sex[pop$sex=="both"] <- "Both"
pop <- pop %>%
  filter(age=="All Ages") %>%
  dplyr::select("location","sex","year","val") %>%
  rename(pop=val)
# 合并人口数据, 生成 mydata
mydata <- left_join(data,pop,
                    by=c("location","sex","year"))


# 斜度指数的可视化 ----------------------------------------------------------------

## 1.绘图数据的准备 -----------------------------------------------------------------
# 计算总人口
a <- mydata %>%
  filter(metric=="Number") %>%
  group_by(year) %>%
  summarise(sum=sum(pop))
pop1990 <- a$sum[1]
pop2019 <- a$sum[2]
# 计算加权次序
rank <- mydata %>%
  mutate(pop_global=ifelse(year==1990,pop1990,pop2019)) %>%
  group_by(year,metric) %>%
  arrange(sdi) %>%
  mutate(cummu=cumsum(pop)) %>% # 计算累积人口
  mutate(half=pop/2) %>% # 计算该国家人口的一半
  mutate(midpoint=cummu-half) %>% # 累积人口减去该国家人口一半即为人口中点
  mutate(weighted_order=midpoint/pop_global) # 人口中点与总人口相比即为改国的相对位置
# 把年份设置为 factor
rank$year <- factor(rank$year)
# 选择数据
temp1 <- rank %>%
  filter(metric=="Rate") %>%
  filter(year==1990)
temp2 <- rank %>%
  filter(metric=="Rate") %>%
  filter(year==2019)
# 建模计算斜度指数
fit1 <- lm(data = temp1,val~weighted_order)
fit2 <- lm(data = temp2,val~weighted_order)
# 查看是否存在异方差（存在异方差）
ncvTest(fit1)
ncvTest(fit2)
# 使用稳健（robust）回归：重复迭代加权
r.huber1 <- rlm(data = temp1,val~weighted_order)
r.huber2 <- rlm(data = temp2,val~weighted_order)
# 获得系数与截距
coef(r.huber1)
coef(r.huber2)
# 计算稳健回归的 95% 可信区间
confint.default(r.huber1)
confint.default(r.huber2)

# 2.绘图 ----------------------------------------------------------------------
color <- c("#6699FF","#990000")
colnames(rank)
p1 <- rank %>%
  filter(metric=="Rate") %>%
  ggplot(aes(x=weighted_order,y=val,fill=year,group=year,color=year))+
  geom_point(aes(color=year,size=pop/1e6),alpha=0.8,shape=21)+
  scale_size_area("Population\n(million)",breaks=c(200,400,600,800,1000,1200))+
  geom_smooth(method = "rlm",size=0.6,alpha=0.1)+
  scale_fill_manual(values = color)+
  scale_color_manual(values = color)+
  # 增加中括号：中括号的位置需要提前确定，根据截距与斜率
  geom_brace(aes(x=c(1.003,1.103),y=c(87548,87548-70136)),
             inherit.data = F,size=0.6,
             rotate=90,color="#6699FF")+
  geom_brace(aes(c(1,1.1),c(40734,40734-13805)),
             inherit.data = F,size=0.6,
             rotate=90,color="#990000")+
  # 增加水平虚线
  geom_segment(x=0.02,xend=0.99,
               y=87548,yend=87548, # 截距的位置
               color="#6699FF",linetype=2,size=0.4,alpha=0.4)+
  geom_segment(x=0.02,xend=0.99,
               y=40734,yend=40734, # 截距的位置
               color="#990000",linetype=2,size=0.4,alpha=0.4)+
  # 增加某些国家的标签: 比如中国与印度
  geom_text(aes(label=ifelse(location=="China"|location=="India",as.character(location),""),
                color=year),
            hjust=0,vjust=1.7,# 避免点和文字重合
            size=3)+
  # 增加斜度指数标签
  annotate("text",label="Slope Index of Inequality",x=1.22,y=40670+6751/2,size=4,angle=90)+
  annotate("text",label="-70136",x=1.15,y=87548-70136/2,size=3.5)+ # 位置
  annotate("text",label="-13805",x=1.15,y=40734-13805/2,size=3.5)+
  scale_x_continuous(limits = c(0,1.22),labels = c("0","0.25","0.50","0.75","1.00",""))+
  xlab("Relative rank by SDI")+
  ylab("Crude DALY rate (per 100,000)")+
  theme_bw()

p1

# 集中指数的可视化 ----------------------------------------------------------------

# 1.绘图数据准备 ------------------------------------------------------------------
a <- mydata %>%
  filter(metric=="Number") %>%
  group_by(year) %>%
  summarise(sum=sum(val))
daly1990 <- a$sum[1]
daly2019 <- a$sum[2]

ci <- rank %>%
  filter(metric=="Number") %>%
  mutate(total_daly=ifelse(year==1990,daly1990,daly2019)) %>%
  group_by(year) %>%
  arrange(sdi) %>%
  mutate(cummu_daly=cumsum(val)) %>% # 计算累积 daly
  mutate(frac_daly=cummu_daly/total_daly) %>% # 计算累积 daly 所占总体的比例
  mutate(frac_population=cummu/pop_global) # 计算累积人口所占总体人口的比例



# 2.绘图 --------------------------------------------------------------------
p2 <- ci %>%
  ggplot(aes(x=frac_population,y=frac_daly,fill=year,color=year,group=year))+
  # 增加 X=0,y=0 两条线段
  geom_segment(x=0,xend=1,
               y=0,yend=0,
               linetype=1,size=1,color="gray")+
  geom_segment(x=1,xend=1,
               y=0,yend=1,
               linetype=1,size=1,color="gray")+
  # 对角线
  geom_segment(x=0,xend=1,
               y=0,yend=1,
               color="#CD853F",linetype=1,size=0.7,alpha=1)+
  # 散点
  geom_point(aes(fill=year,size=pop/1e6),alpha=0.75,shape=21)+
  scale_fill_manual(values = color)+
  scale_size_area("Population\n(million)",breaks=c(200,400,600,800,1000,1200))+
  # 立方样条函数拟合洛伦兹曲线 (设置节点，边界条件)
  geom_smooth(method = "gam", # 这里也可以直接用 geom_line 把点连起来
            formula = y ~ ns(x,
                             knots = c(0.0000000001,0.25,0.5,0.75,0.9999999),# 设置节点为
                             # df = 3, 上面设置了这里就不需要设置
                             Boundary.knots = c(0,1)),
            linetype=1,size=0.1,alpha=0.6,se=T)+
  scale_color_manual(values = color)+
  # 增加两个年份的集中指数
  annotate("text",label="Concentration Index",x=0.75,y=0.35,size=5)+
  annotate("text",label="1990: 0.1?",x=0.75,y=0.3,size=4,color="#6699FF")+
  annotate("text",label="2019: 0.11?",x=0.75,y=0.25,size=4,color="#990000")+
  # 增加某些国家的标签，1990 年
  geom_text(aes(label=ifelse(location=="China"&year=="1990"|location=="India"&year=="1990",
                             as.character(location),"")),
            hjust=-0.6,vjust=0.8,
            size=3)+
  # 增加某些国家的标签，2019 年
  geom_text(aes(label=ifelse(location=="China"&year=="2019"|location=="India"&year=="2019",
                             as.character(location),"")),
            hjust=1.8,vjust=-0.0,
            size=3)+
  # 增加某些国家标签，人口大国
  geom_text(aes(label=ifelse(location%in%a&year=="1990",
                             as.character(location),"")),
            hjust=-0.6,vjust=0.8,
            size=3)+
  geom_text(aes(label=ifelse(location%in%a&year=="2019",
                             as.character(location),"")),
            hjust=1.8,vjust=-0.0,
            size=3)+
  # xy 标签
  xlab("Cumulative fraction of population ranked by SDI")+
  ylab("Cumulative fraction of DALY")+
  theme_bw()
p2



