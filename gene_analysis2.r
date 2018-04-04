#帮助文档非常重要！

#实现数据多项功能检验：

#您需要提供3个数据：data原数据，label哪些列为对照和控制，method您需要采取的检验方法

#1. 对比分析，method=='comparative' || method=='comp' || method==1

#2. 关联分析，method=='association' || method=='asc' || method==2

#3. 相关性分析，method=='correlation' || method=='cor' || method==3

#这是label的定义方法
 #idx.case<-grep('01$',colnames(data),perl=T)
 #label<-rep('control',dim(data)[2])
 #label[idx.case]='case'

gene_analyse<-function(data,label,method){
	alpha=0.05#默认显著性水平为0.05
	data<-as.matrix(data)#防止其它的数据类型对下面计算产生干扰
	
	#1计算是否符合正太分布
	if(dim(data)[1]*dim(data)[2]<3000){
		norm_test<-shapiro.test(data)#适用于小样本3-5000，由于现在的样本数可以说是非常多，所以推荐用ks
	}else{
		norm_test<-ks.test(jitter(data),'pnorm')#此处进行扰动，不会使得
	}
	norm_value<-norm_test$p.value
	
	data1_record<-which(label==levels(as.factor(label))[1])
	data2_record<-which(label==levels(as.factor(label))[2])
	#将其根据label分，认为只有两组，control组与对照组
	
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
		result<-as.matrix(result,ncol=1)
		rownames(result)<-rownames(data)
		
		
	}else if(method=='association' || method=='asc' || method==2){
		result<-c()
		#构造矩阵,对每个gene进行分析
		for(i in 1:dim(data)[1]){
			data_trans<-rep(0,dim(data)[2])#同一个gene在不同样本中，若表达量高于均值则赋为1
			data_trans[which(data[i,]>mean(data[i,]))]=1
			
			x1<-length(which(data_trans[data1_record]==1))
			x2<-length(which(data_trans[data2_record]==1))
			x3<-length(which(data_trans[data1_record]==0))
			x4<-length(which(data_trans[data2_record]==0))
		
			data_association<-matrix(c(x1,x3,x2,x4),nrow=2,dimnames=list(c('high','low'),c('case','control')))
			
			sample_num<-x1+x2+x3+x4
			
			if(sample_num<40){  #样本数小于40或者每个样本的理论频数<1(后面条件自动满足)
				result[i]=fisher.test(data_association)$p.value
			}else{
				result[i]=chisq.test(data_association)$p.value
			}
		}
		result<-p.adjust(result)
		names(result)<-rownames(data)

		
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
