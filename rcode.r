library(Arothron)
library(Morpho)
library(compositions)
library(rgl)
library(alphahull)
library(bezier)
library(Momocs)
library(Rvcg)

all.mesh.set<-function(set,mesh){
  eucl<-dist(set,method="euclidean")
  newP1<-c(0,0,0)
  newP2<-c(eucl[1],0,0)
  newP3<-c(((eucl[1]^2)+(eucl[2]^2)-(eucl[3]^2))/(2*eucl[1]), sqrt((eucl[2]^2)-(((eucl[1]^2)+(eucl[2]^2)- (eucl[3]^2))/(2*eucl[1]))^2), 0)
  tar<-rbind(newP1,newP2,newP3)
  rot_mesh<-rotmesh.onto(mesh,as.matrix(set),as.matrix(tar))
  return(rot_mesh)
}
order_matrix=function(matrix,v1,v2){
  Vector=sqrt(c(v1-matrix[,1])^2+c(v2-matrix[,2])^2)
  pos=c(which(Vector==min(Vector)))
  for(i in 1:(nrow(matrix)-1)){
    vector=sqrt(c(matrix[pos[i],1]-matrix[,1])^2+c(matrix[pos[i],2]-matrix[,2])^2)
    vector_junk=vector[-c(which(vector==min(vector)),pos)]
    pos=c(pos,which(vector==sort(vector_junk)[1]))
  }
  matrix_out=matrix[pos,]
  return(matrix_out)
}
sort_points <- function(df,clockwise = TRUE) {
  # Get centre (-oid) point of points
  x_centre <- mean(df[, 1])
  y_centre <- mean(df[, 2])
  
  # Calculate deltas
  df$x_delta <- df[, 1] - x_centre
  df$y_delta <- df[, 2] - y_centre
  
  # Resolve angle, in radians
  df$angle <- atan2(df$y_delta, df$x_delta)
  # d$angle_degrees <- d$angle * 180 / pi
  
  # Arrange by angle
  if (clockwise) {
    
    df <- df[order(df$angle, decreasing = TRUE), ]
    
  } else {
    
    df <- df[order(df$angle, decreasing = FALSE), ]
    
  }
  
  # Drop intermediate variables
  df[, c("x_delta", "y_delta", "angle")] <- NULL
  
  # Return
  df
  
}

##import coordinates along the midsagittal plane
#set_mds<-read.amira.set("mds.txt","auto")[,,1]

##import femure bone
#all_sur<-file2mesh("femur_left_ori.ply")
rot_all_sur<-all.mesh.set(set_mds,all_sur)

bbox<-meshcube(rot_all_sur$mesh)
diff_length<-abs(bbox[1,1])
rot_all_sur$mesh$vb[1:3,]=(rot_all_sur$mesh$vb[1:3,]+c(diff_length,0,0))
#mesh2ply(rot_all_sur$mesh,"mesh_all")

ca_lse<-ext.int.mesh(rot_all_sur$mesh, views=30, param1=3, default=TRUE,
                     import_pov = NULL)
out_inn<-out.inn.mesh(ca_lse,rot_all_sur$mesh)                          

inner_sur<-out_inn$invisible
outer_sur<-out_inn$visible

wire3d(vcgQEdecim(outer_sur,10000),col="white")
shade3d(inner_sur,col="violet")

bio_length<-as.matrix(dist(meshcube(rot_all_sur$mesh)))[1,2]
sect_poi<-seq(bio_length*0.2,bio_length*0.8,length=61)
out_coo_3D<-array(NA,dim=c(73,3,length(sect_poi)))
inn_coo_3D<-array(NA,dim=c(73,3,length(sect_poi)))
out_coo_2D<-array(NA,dim=c(21,2,length(sect_poi)))
inn_coo_2D<-array(NA,dim=c(21,2,length(sect_poi)))

dist_sections<-array(NA,dim=c(73,1,length(sect_poi)))
for(m in 1:length(sect_poi)){
  p1<-c(sect_poi[m],0,0)
  p2<-c(sect_poi[m],20.5,0.5)
  p3<-c(sect_poi[m],-20.5,1.5)
  inter_out<-NULL
  inter_inn<-NULL
  inter_out <- meshPlaneIntersect(outer_sur,p1,p2,p3)
  inter_inn <- meshPlaneIntersect(inner_sur,p1,p2,p3)
  
  spheres3d(inter_out)
  spheres3d(inter_inn,col=2)
  
  coord_out<-ashape(inter_out[,2],inter_out[,3],alpha=30)
  coord_inn<-ashape(inter_inn[,2],inter_inn[,3],alpha=30)
  
  mat_out<-cbind(inter_out[coord_out$alpha.extremes,2],inter_out[coord_out$alpha.extremes,3])
  mat_inn<-cbind(inter_inn[coord_inn$alpha.extremes,2],inter_inn[coord_inn$alpha.extremes,3])
  
  # plot(mat_out,asp=1)
  # points(mat_inn)
  ordered_out_temp<-sort_points(as.data.frame(mat_out))
  # 
  # ordered_out_temp<-order_matrix(mat_out,mat_out[which.max(mat_out[,2]),1],
  #                                mat_out[which.max(mat_out[,2]),2])
  # ordered_out_temp<-order_matrix(mat_out,apply(mat_out,2,mean)[1],
  #                                apply(mat_out,2,mean)[2])[-1,]
  ordered_out<-rbind(ordered_out_temp,ordered_out_temp[1,])
  # plot(ordered_out,asp=1)
  # text(ordered_out,labels=c(1:dim(ordered_out)[1]))
  ordered_inn_temp<-sort_points(as.data.frame(mat_inn))
  
  # ordered_inn_temp<-order_matrix(mat_inn,mat_out[which.max(mat_inn[,2]),1],
  #                                mat_out[which.max(mat_inn[,2]),2])
  ordered_inn<-rbind(ordered_inn_temp,ordered_inn_temp[1,])
  
  
  out_bez<-pointsOnBezier(ordered_out,method = "max_dist",
                          max.dist = 0.1)$points
  inn_bez<-pointsOnBezier(ordered_inn,method = "max_dist",
                          max.dist = 0.1)$points
  ev_out<-equidistantCurve(out_bez,n=22,iterations=10,increment=0)
  ev_inn<-equidistantCurve(inn_bez,n=22,iterations=10,increment=0)
  
  out_coo_2D[,,m]<-ev_out[-1,]
  inn_coo_2D[,,m]<-ev_inn[-1,]
  
  # out_bez<-out_bez[round(seq(1,dim(out_bez)[1],length=100)),]
  # inn_bez<-inn_bez[round(seq(1,dim(inn_bez)[1],length=100)),]
  
  # plot(out_bez,asp=1,type="l")
  # points(inn_bez,type="l",col="red")
  
  coo_i <- inn_bez[-1,]
  # coo_plot(coo_i)
  sapply(seq(0, 2*pi, pi/36),
         function(x) coo_i %>% coo_intersect_angle(x)) -> ids_i
  # coo_i[ids_i, ] %>% points(col="blue")
  
  coo_o <- out_bez[-1,]
  # coo_plot(coo_o)
  ids_o<-NULL
  sapply(seq(0, 2*pi, pi/36),
         function(x) coo_o %>% coo_intersect_angle(x)) -> ids_o
  # coo_o[ids_o, ] %>% points(col="blue")
  
  # aaa<-NULL
  # for(x in seq(0, 2*pi, pi/36)){
  # aaa<-c(aaa,coo_intersect_angle(coo_o[2:20,],round(x,1)))
  # }
  
  
  
  # plot(coo_o[ids_o, ],asp=1,col="red",type="p")
  # points(coo_i[ids_i, ],asp=1,col="black",
  #        type="p")
  
  dist_section<-sqrt(rowSums((coo_o[ids_o, ]-coo_i[ids_i, ])^2))
  cols<-rainbow(73)
  dist_sections[,,m]<-dist_section
  ramp <- colorRampPalette(c("green","orange","red"))(100)
  ord_cols<-order(dist_section)
  
  intervals<-seq(min(dist_section),max(dist_section),length=100)
  pos_cols<-NULL
  for(i in 1:length(dist_section)){
    pos_cols[i]<-which.min(abs(dist_section[i]-intervals))
  }
  
  
  plot(coo_o[ids_o, ],asp=1,col=ramp[pos_cols],
       type="p",pch=19,cex=3,main=paste("section at ",19+m,"%",sep=""),
       xaxt='n',yaxt='n',xlab="",ylab="")
  points(coo_i[ids_i, ],asp=1,col="black",
         type="l")
  
  A<-coo_o[ids_o, ]
  B<-coo_i[ids_i, ]
  for(i in 1:length(dist_section)){
    points(rbind(A[i,],B[i,]),type="l",col="violet")
  }
  
  out_coo_3D[,,m]<-cbind(sect_poi[m],coo_o[ids_o, ])
  inn_coo_3D[,,m]<-cbind(sect_poi[m],coo_i[ids_i, ])
}

aaaa<-vcgClostKD(out_coo_3D[,,1],outer_sur)

plot3D(t(aaaa$vb)[,-4])
matricione<-out_coo_3D[,,1]
for(i in 2:dim(out_coo_3D)[3]){
  matricione<-rbind(matricione,out_coo_3D[,,i])  
}

plot3D(matricione)
spess_values<-as.vector(unlist(dist_sections))
dim(dist_sections)
cols<-colorRampPalette(c("green","orange","red"))(100)
intervals<-seq(min(spess_values),max(spess_values),length=100)
pos_cols<-NULL
for(i in 1:length(spess_values)){
  pos_cols[i]<-which.min(abs(spess_values[i]-intervals))
}

open3d()
spheres3d(matricione,col=cols[pos_cols])
wire3d(vcgQEdecim(outer_sur,10000),col="white")

# matri<-t(cbind(matricione,1))
# surm<-list("vb"=matri,"it"=NULL)
# class(surm)="mesh3d"
# aaa<-vcgClost(outer_sur,surm)
# str(aaa)
# dim(matricione)
# shade3d(surm)

clostInd <- mcNNindex(matricione,t(outer_sur$vb)[,1:3]
                      ,k=20)
distInd<-clostInd
for(i in 1:ncol(clostInd)){
  distInd[,i]<-sqrt(rowSums((t(outer_sur$vb)[,1:3]-matricione[clostInd[,1],])^2))
}

colsInd<-NULL
for(i in 1:dim(clostInd)[1]){
  colsInd[i]<-round(weighted.mean(pos_cols[clostInd[i,]], distInd[i,]/sum(distInd[i,])))
}

white_col<-c(which(t(outer_sur$vb)[,1]<bio_length*0.2),
which(t(outer_sur$vb)[,1]>bio_length*0.8))

cols<-c(cols,"#FFFFFF")
colsInd[white_col]<-101
mesh2ply(outer_sur,filename = "prova",col=cols[colsInd])
aaa<-file2mesh("prova.ply",readcol = TRUE)
shade3d(aaa)
spheres3d(lds,col=1,radius=2)

open3d()
spheres3d(matricione,col=cols[pos_cols])

library(grDevices)
library(graphics)




filled.contour(matrix(spess_values,nrow=73,ncol=61),
               color.palette = colorRampPalette(c("green","orange","red")), 
               nlevels=10,asp=1)

##import set to be plotted on the 2D map
#path_lines<-read.path.amira("path_aspra.am")
str(path_lines)
poi3D<-pointsOnBezier(unique(path_lines[seq(1,dim(path_lines)[1],length.out = 100),]),11)
lds<-poi3D$points

z1 <- apply(out_coo_3D, 2L, c)
landpos<-mcNNindex(z1,lds,k=1)

arr1pos<-NULL
arr2pos<-NULL
for(i in 1:length(landpos)){
arr1pos[i]<-which(out_coo_3D[,1,][1,]==z1[landpos[i],1])
arr2pos[i]<-which(out_coo_3D[,2,arr1pos[i]]==z1[landpos[i],2])
}

filled.contour(matrix(spess_values,nrow=73,ncol=61),
color.palette = colorRampPalette(c("green","orange","red")), 
nlevels=10,asp=1,plot.axes={points(cbind(arr2pos/73,arr1pos/61),pch=19)})

  
