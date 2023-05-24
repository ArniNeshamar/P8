
library(ggplot2)


DistMat <- function(locations, method = "euclidean"){
  # locations -- [ID, X, Y] Node id and X, Y co-ordinates
  # method = "euclidean"
  if(method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")){
    X = locations[,2]%*%t(rep(1,length(locations[,1])));  Y = locations[,3]%*%t(rep(1,length(locations[,2])))
    DMat <- sqrt((X-t(X))^2+(Y-t(Y))^2)
  }
  DMat <- as.matrix(DMat)
  row.names(DMat) <- locations[, 1]
  colnames(DMat) <- locations[, 1]
  return(DMat)
}

Total_Cost <- function(result, DMat){
  cost <- 0
  for(i in 1:length(result)){
    for(j in 1:(length(result[[i]]) - 1)){
      cost <- cost + DMat[result[[i]][j], result[[i]][j+1]]
    }
  }
  return(cost)
}

SavingMat <- function(DMat, Routes = NULL, Weights = NULL, Wei = FALSE, depot = 1){
  SMat <- matrix(0, nrow = nrow(DMat) - 1, ncol = nrow(DMat) - 1)
  row.names(SMat) <- row.names(DMat)[-1]
  colnames(SMat) <- row.names(DMat)[-1]
  dist = c()
  for(k in 1:length(Routes)){
    d = DMat[1,Routes[[k]][length(Routes[[k]])]]
    if(length(Routes[k])>1){
      for(i in 2:length(Routes[[k]])){
        d = d + DMat[Routes[[k]][i-1],Routes[[k]][i]]
      }
    }
    dist = c(dist,d)
  }
  for(i in seq_len(nrow(SMat))){
    w1 = Inf
    for(k in 1:length(Routes)){
      if(Routes[[k]][length(Routes[[k]])] == (i+1)){ w1 = Weights[[k]] }
    }
    for(j in seq_len(ncol(SMat))){
      w2 = Inf;  j1 = 1
      for(k in 1:length(Routes)){
        if(Routes[[k]][1] == (j+1)){ w2 = Weights[[k]];  j1 = k }
      }
      if(i != j){
        if(!Wei){ SMat[i, j] = DMat[depot, i+1] + DMat[depot, j+1] - DMat[i+1, j+1] }
        else if(w1 == Inf){ SMat[i,j] = 0 }
        else if(w2 == Inf){ SMat[i,j] = 0 }
        else{
          SMat[i, j] = (1+0.01*w1)*DMat[depot, i+1] + DMat[depot, j+1] - (1+0.01*w1)*DMat[i+1, j+1] - (0.01*w1)*dist[j1]
        }
      }
    }
  }
  return(SMat)
}

Sorted_Edges <- function(DMat, Routes = NULL, Weights = NULL, Wei = FALSE){
  SMat <- SavingMat(DMat, Routes = Routes, Weights = Weights, Wei = Wei, depot = 1)
  result <- as.data.frame.table(SMat)
  colnames(result) <- c("i", "j", "Saving")
  result$i <- as.integer(as.character(result$i))
  result$j <- as.integer(as.character(result$j))
  result <- result[result$Saving != 0, ]
  result <- result[order(result$Saving, decreasing = TRUE), ]
  result$Saving <- as.numeric(result$Saving)
  return(result)
}

Sorted_Edges(DMat)



CW = function(Locations,DMat,Demand,Wei=FALSE,Cluster = rep(1,length(Locations[,1])-1),title = "Plot of Greedy Routes", capacity = 100,plot=TRUE){
  n = length(DMat[,1])-1
  routes = list()
  weights = list()
  for(i in 1:n){
    routes = append(routes,list(c(i+1)))
    weights = append(weights,list(c(Demand[i+1])))
  }
  W = TRUE
  while(W){
    SE = Sorted_Edges(DMat,Routes=routes,Weights=weights, Wei = Wei)
    for(i in 1:length(SE[,1])){
      N = c(SE[i,1],SE[i,2])
      Avail = c(0,0)
      for(j in 1:length(routes)){
        if(routes[[j]][length(routes[[j]])] == N[1]){ Avail[1] = j }
        else if(routes[[j]][1] == N[2]){ Avail[2] = j }
      }
      if(prod(Avail) != 0){
        if(weights[[Avail[1]]]+weights[[Avail[2]]] >= capacity){ next }
        R1 = routes[[Avail[1]]];  R2 = routes[[Avail[2]]]
        routes = append(routes[-Avail],list(c(R1,R2)))
        w = weights[[Avail[1]]] + weights[[Avail[2]]]
        weights = append(weights[-Avail],w)
        if(Wei){ break }
      }
    }
    if(i == length(SE[,1])){ W = FALSE }
  }
  for(i in 1:length(routes)){
    routes[[i]] = c(1,routes[[i]],1)
  }
  if(plot){
    g <- ggplot(Locations[unlist(routes), ], aes_string(x = names(Locations)[2], y = names(Locations)[3])) + geom_path(lineend = "round", linetype = 2, show.legend = TRUE) + labs(title = title) 
    g <- g + annotate("text", x = Locations[, 2], y = Locations[, 3], label = 1:(n+1), col=c(1,1+Cluster))
    plot(g)
  }
  #print(routes)
  #print(cumsum(Demand[routes[[1]]][c(-1,-51)]))
  #print(cumsum(Demand[rev(routes[[1]][c(-1,-51)])]))
  return(routes)
}

GenerateLocations = function(n){
  loc = matrix(c(1:n,round(100*runif(2*n))),n,3)
  rownames(loc) = 1:n;  colnames(loc) = c("ID","X","Y")
  loc = as.data.frame(loc)
  return(loc)
}

GenerateDemand = function(n){
  loc = matrix(c(1:n,round(2*rgamma(n,5))),n,2)
  rownames(loc) = 1:n;  colnames(loc) = c("ID","Demand")
  loc = as.data.frame(loc)
  return(loc)
}

locations = GenerateLocations(50)
demand = GenerateDemand(50)
t(demand)
DMat = DistMat(locations)

Total_Cost(CW(locations,DMat,demand$Demand,Cluster = Cluster,title = "Standard CW",capacity = 1000),DMat)
Total_Cost(CW(locations,DMat,demand$Demand,Cluster = Cluster,title = "Standard CW with capacity"),DMat)
Total_Cost(CW(locations,DMat,demand$Demand,TRUE,Cluster = Cluster, title = "Standard CW with Weight",capacity=1000),DMat)
Total_Cost(CW(locations,DMat,demand$Demand,TRUE,Cluster = Cluster, title = "Standard CW with capacity and Weight"),DMat)
demand$Demand
list(1,2,3)[-1]

CW.mean = 0
LP.mean = 0
CapDist = c()

for(i in 87:100){
  L = GenerateLocations(50)
  D = GenerateDemand(50)
  DM = DistMat(L)
  C = kmeans(L[-1,2:3],8)[[1]]
  #CW.mean = CW.mean + Total_Cost(CW(L,DM,D$Demand, Cluster = C,title = "Standard CW with capacity",plot=F),DM)
  R = ClusterSingleRoute(L,C,title="PP",plot=F)
  for(j in 1:length(R)){
    CapDist = c(CapDist,sum(D$Demand[R[[j]][2:(length(R[[j]])-1)]]))
  }
  LP.mean = LP.mean + Total_Cost(R,DM)
}
CW.mean/100
LP.mean/100
hist(CapDist,20,xlab="Garbage Load",main="Histogram of Garbage Load in Cluster LP")
mean(CapDist>100)

{
  L[1,]$X = 50;  L[1,]$Y = 50
  DM = DistMat(L);  D$Demand = 0
  TotalDist = 0;  TrashHist = c()
  td = 0;  tn = 0;  dn = 0
  deadlines = rep(10,50)
  for(i in 1:100){
    D$Demand = D$Demand + rgamma(50,2,1)
    D[1,]$Demand = 100
    TrashHist = c(TrashHist,D[-1,]$Demand)
    if(max(D[-1,]$Demand)>10){ if(td < i){ td = i } }
    if(td == i){
      dn = dn + 1
      D$Demand = round(D$Demand)
      Ind = (1:50)[(D$Demand > 10)|(deadlines == 0)]
      if(length(Ind) < 3){ deadlines = deadlines - 1;  next }
      print(i)
      print(Ind)
      DM.temp = DM[Ind,Ind]
      colnames(DM.temp) = 1:length(Ind)
      rownames(DM.temp) = 1:length(Ind)
      Route = CW(L[Ind,],DM.temp,D[Ind,]$Demand,title = paste0("Standard CW with capacity Simulation, day ",i),plot=F)
      Route = OnTheWay(Route,Ind,L);  Ind = sort(unique(unlist(Route)))
      print(Ind);  print(Route)
      D[Ind,]$Demand = 0
      TotalDist = TotalDist + Total_Cost(Route,DMat)
      tn = tn + length(Ind)
      deadlines[Ind] = 10
    }
    deadlines = deadlines - 1
  }
}

hist(TrashHist,main="Histogram of Daily Amount of Trash in Bins, Pickup at 10", xlab = "Amount of Trash")
mean(TrashHist>12)
TotalDist
tn
dn
locations

OnTheWay = function(Routes,Ind,Locations,title="",threshhold=1){
  NewRoutes = Routes
  for(r in 1:length(Routes)){
    Route = Routes[[r]]
    Route = Ind[Route];  NewRoute = Route
    n = length(Locations[,1]);  l = length(Route);  I = 0
    for(j in 2:l){
      m = Inf; ind = 0
      for(i in 1:n){
        if(i %in% Route){ next }
        D = (locations$X[i]-Locations$X[Route[j]])^2 + (Locations$Y[i]-Locations$Y[Route[j]])^2
        D = D + (Locations$X[i]-Locations$X[Route[j-1]])^2 + (Locations$Y[i]-Locations$Y[Route[j-1]])^2
        D = D - (Locations$X[Route[j]]-Locations$X[Route[j-1]])^2 + (Locations$Y[Route[j]]-Locations$Y[Route[j-1]])^2
        if(D < m){ ind = i;  m = D }
      }
      if(m < threshhold){ NewRoute = c(NewRoute[1:(j+I-1)],ind,NewRoute[(j+I):(l+I)]);  I = I + 1 }
    }
    NewRoutes[[r]] = NewRoute
  }
  if(F){
    g <- ggplot(Locations[unlist(Routes), ], aes_string(x = names(locations)[2], y = names(Locations)[3])) + geom_path(lineend = "round", linetype = 2, show.legend = TRUE) + labs(title = title) 
    g <- g + annotate("text", x = Locations[, 2], y = Locations[, 3], label = 1:n, col=1)
    plot(g)
    g <- ggplot(Locations[unlist(NewRoutes), ], aes_string(x = names(locations)[2], y = names(Locations)[3])) + geom_path(lineend = "round", linetype = 2, show.legend = TRUE) + labs(title = title) 
    g <- g + annotate("text", x = Locations[, 2], y = Locations[, 3], label = 1:n, col=1)
    plot(g)
  }
  return(NewRoutes)
}

OnTheWay(Route,Ind,L)

points(locations[Ind[Route[[1]]],2:3])
plot(locations[OnTheWay(Route[[1]],Ind,locations),2:3],col="red")
Route


