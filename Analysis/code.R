setwd("~/Desktop/Bayes")

library(doParallel)
library(foreach)
library(raster)
library(sp)
library(rgdal)
library(ggmap)
library(plotly)
library(foreign)
library(ggplot2)

############################Extract satellite raw image, based on Star Ying's code
######Set directory path
imagery = "imagery"

######Obtain a list of TIF files, load in the first file in list
tifs = list.files(imagery,pattern = "\\.tif")
rast <- raster(paste(imagery,"/",tifs[1],sep=""))

######Specify WGS84 as the projection of the raster file
wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
projection(rast) <- CRS(wgs84)

######Draw down county Shapefile
shape_direct <- function(url, shp) {
  library(rgdal)
  temp = tempfile()
  download.file(url, temp) ##download the URL taret to the temp file
  unzip(temp,exdir=getwd()) ##unzip that file
  return(readOGR(paste(shp,".shp",sep=""),shp))
}

###### download county shapefile
msa <- shape_direct(url="http://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_county_20m.zip", 
                    shp= "cb_2014_us_county_20m")
projection(msa) <- CRS(wgs84)


###read in county id
county<-read.csv("county.csv",header=T) ##a table with county and state id
county$check<-county$FIPS.Code/100
county<-county[county$FIPS.Code%%100!=0,]
county$name<-paste(county$Area.name,county$State,sep=",")

###############################Extract radiance using county shape file
masq <- function(shp,rast,i){
  
  #Extract one polygon based on index value i
  polygon <- shp[i,] #extract one polygon
  extent <- extent(polygon) #extract the polygon extent 
  
  #Raster extract
  outer <- crop(rast, extent) #extract raster by polygon extent
  inner <- mask(outer,polygon) #keeps values from raster extract that are within polygon
  data <- as.matrix(inner)
  data[is.na(data)] = 0 ##replace NA with no light = 0
  return(data)
}

######run 
mat<-list()
for (i in 1:dim(msa)[1]){
  mat[[i]] <- masq(msa,rast,i)
}

############################## project everybody to median image with size 103 131
ref<-raster(mat[[1212]])

re_result<-lapply(1:3220, function (i) resample(raster(mat[[i]]),ref,method="bilinear"))

######convert raster to vector
re_result2<-lapply(1:3220, function (i) as.vector(re_result[[i]]))
normalize<-matrix(unlist(re_result2), nrow = 3220)

######remove columns with all 0 since later, we will do PCA which requires centering and scaling
new<-normalize[ , !apply(normalize==0,2,all)] 
image.pca <- prcomp(new,center = TRUE,scale. = TRUE) ## PCA for dim reduction

newT<-t(normalize) ##to compute eigen county, common county patterns
dim(newT)
image.pca4 <- prcomp(newT,center = T,scale = T) ##PCA for eigen county

scores<-data.frame(image.pca$x) ##use this for dim reduction, projected data on principal components
scores4<-data.frame(image.pca4$x)##use this for eigen city

######project data onto first two principles
ggplot(data=scores2,aes(x=PC1,y=PC2,label=rownames(scores2)))+
  geom_hline(yintercept=0,colour="gray65")+
  geom_vline(xintercept=0,colour="gray65")+
  geom_text(colour="tomato",alpha=0.2,size=4)+
  ggtitle("PCA Plot of First Two Principal Components")

plot(image.pca)

######use first 8 components based on the above eigen decay plot

princ_scores<-scores4[1:8]
######select first 4 eigen counties
pc1<-data.matrix(princ_scores[1])
pc2<-data.matrix(princ_scores[2])
pc3<-data.matrix(princ_scores[3])
pc4<-data.matrix(princ_scores[4])
pc1<-matrix(pc1,nrow=103)
pc2<-matrix(pc2,nrow=103)
pc3<-matrix(pc3,nrow=103)
pc4<-matrix(pc4,nrow=103)

######plot eigen county
par(mai=c(0.2,0.2,0.2,0.2),mfrow = c(1,1),bg='white', bty='n')
sampled <- as.vector(pc1)
clusters <- 30
clust <- kmeans(sampled,clusters,iter.max = 1000)$cluster
combined <- as.data.frame(cbind(sampled,clust))
brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])
plot(raster(pc1), breaks=brk,col=colorRampPalette(c("#001a4d","#0066FF", "blue"))(clusters), 
     legend=F,yaxt='n',xaxt='n',frame = F, asp=1.5)


######Perform kmeans on reduced dim data
wss<-(nrow(princ_scores)-1)*sum(apply(princ_scores,2,var))
for (i in 1:20) wss[i]<-sum(kmeans(princ_scores,centers=i,iter.max=1000)$withinss)
plot(1:20,wss,type="b",xlab="Number of Clusters",ylab="Within Groups Sum of Squares")


######Use K=9 in clustering
princ_scores$cluster<-kmeans(princ_scores,centers=9,iter.max=1000)$cluster



###############compare original image and projected image
par(mai=c(0.2,0,0,0),mfrow = c(2,4),bg='white', bty='n')
for (i in 1:4){
  sampled <- as.vector(mat[[i]])
  clusters <- 25
  clust <- kmeans(sampled,clusters,iter.max = 1000)$cluster
  combined <- as.data.frame(cbind(sampled,clust))
  brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])
  plot(raster(mat[[i]]), breaks=brk,col=colorRampPalette(c("#001a4d","#0066FF", "yellow"))(clusters), 
       legend=F,yaxt='n',xaxt='n',frame = F, asp=NA)
}

for (j in 1:4){
  
  sampled <- as.vector(re_result[[j]])
  clusters <- 25
  clust <- kmeans(sampled,clusters,iter.max = 1000)$cluster
  combined <- as.data.frame(cbind(sampled,clust))
  brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])
  plot(re_result[[j]], breaks=brk,col=colorRampPalette(c("#001a4d","#0066FF", "yellow"))(clusters), 
       legend=F,yaxt='n',xaxt='n',frame = F, asp=NA,title="ddddd")
}


########stability based K-means clustering, decide number of clusters
get_stable<-function(data,k,m){##k is number of clusters, m is sample proportion
  
  size<-m*nrow(data)
  data$id<-1:nrow(data)
  
  ##subsample
  sample1_id<-sample(1:nrow(data), size, replace = FALSE)
  sample2_id<-sample(1:nrow(data), size, replace = FALSE)
  
  sub_1<-data[sample1_id, ]
  sub_2<-data[sample2_id, ]
  
  ##k-means
  result_1<-data.frame(kmeans(sub_1[2:9], centers = k,algorithm= "Lloyd",iter.max = 1000)$cluster) ##try to alleviate the impact of random starting value
  result_2<-data.frame(kmeans(sub_2[2:9], centers = k,algorithm= "Lloyd",iter.max = 1000)$cluster)
  
  ##cluster result
  sub_1_clust<-data.frame(cbind(sub_1$id,result_1))
  names(sub_1_clust)[1] <- "id"
  names(sub_1_clust)[2] <- "cluster"
  sub_2_clust<-data.frame(cbind(sub_2$id,result_2))
  names(sub_2_clust)[1] <- "id"
  names(sub_2_clust)[2] <- "cluster"
  
  ##find intersected obs and its clustering membership, ensure 1-1 match
  intersect<-merge(sub_1_clust,sub_2_clust,by=c("id"))
  clust_1<-as.integer(as.character(intersect$cluster.x))
  clust_2<-as.integer(as.character(intersect$cluster.y))
  
  ##construct C_ij
  inter.dim <- dim(intersect)[1]
  C_1 <- matrix(clust_1, nr = inter.dim, nc = inter.dim) == matrix(clust_1, nr = inter.dim, nc = inter.dim, byrow = TRUE)
  C_2 <- matrix(clust_2, nr = inter.dim, nc = inter.dim) == matrix(clust_2, nr = inter.dim, nc = inter.dim, byrow = TRUE)
  diag(C_1) <- 0
  diag(C_2) <- 0
  
  ##compute similarity measure
  jaccard <- sum(C_1 * C_2)/(sum(C_1) + sum(C_2) - sum(C_1 * C_2))
  matching<- (sum(C_1 * C_2)+sum((1-C_1) * (1-C_2)))/(sum(C_1 * C_2)+sum((1-C_1) * (1-C_2))+sum((1-C_1)*C_2)+sum((1-C_2)*C_1))
  corr<-sum(C_1 * C_2)/sqrt(sum(C_1)*sum(C_2))
  return(c(jaccard,matching,corr))
}


stable_result_jac<-matrix(0,nrow=14,ncol=100)
stable_result_mat<-matrix(0,nrow=14,ncol=100)
stable_result_cor<-matrix(0,nrow=14,ncol=100)


##################run stability-based K-means
for (k in 2:15){
  for (i in 1:100){
    dd<-get_stable(withname,k=k,m=0.6)
    stable_result_jac[k-1,i]<-dd[1]
    stable_result_mat[k-1,i]<-dd[2]
    stable_result_cor[k-1,i]<-dd[3]
  }
}

jac<-data.frame(t(stable_result_jac))
mat2<-data.frame(t(stable_result_mat))
cor<-data.frame(t(stable_result_cor))

######store results
write.csv(jac,"jac_result.csv")
write.csv(mat2,"mat_result.csv")
write.csv(cor,"cor_result.csv")


#######plot results of matching metric
mat2<-read.csv("mat_result.csv",header=T)
cor<-mat2
df <- data.frame(x = c(cor$X1,cor$X2,cor$X3,cor$X4,cor$X5,cor$X6,cor$X7,cor$X8,cor$X9,cor$X10,cor$X11,cor$X12,cor$X13,cor$X14),ggg = factor(rep(1:14, c(rep(100,14)))))
df <- df[order(df$x), ]
df$ecdf <- ave(df$x, df$ggg, FUN=function(x) seq_along(x)/length(x))

ggplot(df, aes(x, ecdf, colour = ggg)) + geom_line() + 
  scale_colour_hue(name="", labels=c('k=2','k=3','k=4','k=5','k=6','k=7','k=8','k=9','k=10','k=11','k=12','k=13','k=14',"k=15"))+
  xlab("Similarity")+ylab("Cumulative Density")+
  theme(axis.title.y=element_text(family="serif",size=12))+
  theme(axis.title.x=element_text(family="serif",size=12))+
  scale_x_continuous(limits=c(0.4,1))

