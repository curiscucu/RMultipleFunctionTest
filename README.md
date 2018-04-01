# RMultipleFunctionTest
实现数据多项功能检验： 选择数据分析类型
1. 对比分析，判断数据是否服从正态分布： 是，则采用t检验 否，则采用wilcoxon秩和检验 
2. 关联分析，判断数据大小，及四格列联表数值 若总体数据量&lt;40，或单格数据小于1，fisher检验 否则，卡方检验 
3. 相关性分析，判断数据是否服从正态分布： 是，Pearson相关分析 否，Spearman相关分析 
对比分析需要对结果p值进行多重检验校正 
输出结果

对于TCGA 癌症样本 和 正常样本作为conrol 和 Test ,对于这个在癌症中的区分方法有两种，都放在了github上
后续会把这个function包装成shiny app
