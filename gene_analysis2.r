#帮助文档非常重要！

#实现数据多项功能检验：

#您需要提供3个数据：data原数据，label哪些列为对照和控制，method您需要采取的检验方法

#1. 对比分析，method=='comparative' || method=='comp' || method==1

#2. 关联分析，method=='association' || method=='asc' || method==2

#3. 相关性分析，method=='correlation' || method=='cor' || method==3

 

gene_analyse<-function(data,method){
	alpha=0.05#默认显著性水平为0.05 如果p.value大于0.05 则服从正态分布
	data<-as.matrix(data)#防止其它的数据类型对下面计算产生干扰
	
	#1计算是否符合正太分布
	if(dim(data)[1]*dim(data)[2]<3000){
		norm_test<-shapiro.test(data)#适用于小样本3-5000，由于现在的样本数可以说是非常多，所以推荐用ks，
	}else{
		norm_test<-ks.test(jitter(data),data,pnorm)#因为数据不能有重复值，但是数据里确实有重复值，要用jitter()做一个小小的扰动
	}
	norm_value<-norm_test$p.value
	
	#分为control test组
	#这里由于是TCGA数据，分类根据barcode TCGA.02.0321.01 最后这个数字 	Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29.

	group<-c()
	
	for(i in 1:dim(data)[2])
	{
		Split<-unlist(strsplit(colnames(data)[i],"[.]")) 
		group[i]=Split[4]
		
	}
	newData<-rbind(data,group)#newData[12043,]为标签#此处可改进的 

	num<-dim(newData)[1]

	
	data1_record<-which("11"==newData[num,])
	data2_record<-which("01"==newData[num,])


	result<-c()
	
	#2各种检验
	if(method=='comparative' || method=='comp' || method==1){
		if(norm_value>alpha){
			for(i in 1:dim(data)[1]){
				result[i]=t.test(data[i,data1_record],data[i,data2_record])$p.value
			}
		}
		else{
			for(i in 1:dim(data)[1]){
				result[i]=wilcox.test(data[i,data1_record],data[i,data2_record])$p.value
			}
		}
		
		
		result<-p.adjust(result)
	}else if(method=='association' || method=='asc' || method==2){
		data_ass<-apply(data,2,mean)#将样本转化成0,1矩阵的形式，认为每个样本若gene表达的均值高于总体均值则为高表达，低于则低表达
	
		data_association<-matrix(c(length(which(data_ass[data1_record]>mean(data)))
								,length(which(data_ass[data2_record]>mean(data)))
								,length(which(data_ass[data1_record]<mean(data)))
								,length(which(data_ass[data2_record]<mean(data)))),
								nrow=2,
								dimnames=list(c('high','low'),c('case','control'))
								)
								
		sample_num<-dim(data)[1]*dim(data)[2]
		if(sample_num<40 || length(which(data>1))<0){  #样本数小于40或者每个样本的理论频数<1
			result=fisher.test(data_association)
		}
		else{
			result=chisq.test(data_association)
		}
		
	
	}else if(method=='correlation' || method=='cor' || method==3){
		if(norm_value>alpha){
			 result<-cor(t(data),method='pearson')
		}
		else{
			 result<-cor(t(data),method='spearman')
		}
		
	}else{
		stop('What you input has questions!');
	}
	
	
	
	#3
	return (result)
}