# Title     : BinsegRcpp Neuroblastoma Analysis
# Objective : Analyse the Neuroblastoma dataset in order to find the changepoints using the BinsegRcpp package.
# Created by: diego.urgell
# Created on: 13/02/21

packages <- c("neuroblastoma", "ggplot2", "binsegRcpp", "data.table", "devtools", "changepoint", "gfpop")
found <- packages %in% installed.packages()
install.packages(packages[!found])

for(pack in packages) library(pack, character.only=TRUE)


options(warn=-1)

data(neuroblastoma)

nb.profiles = neuroblastoma[['profiles']]
filter <- nb.profiles[['profile.id']] == 4 & nb.profiles[['chromosome']] == 2
nb.profiles <- data.frame(nb.profiles[filter,])

plot(nb.profiles[['position']], nb.profiles[['logratio']], type="l", col="blue",
     xlab="Position in base pairs", ylab="Normalized logratio")
grid()
mtext(line=2, cex=1.2, adj=0.5, expression(bold("Neuroblastoma: Probe logratio vs Position")))
mtext(line=1, cex=1, adj=0.5, "profile.id = 4, chromosome = 2")


models.binsegRcpp <- binseg_normal(nb.profiles[['logratio']], 5)

models.binsegRcpp.coef <- coef(models.binsegRcpp)

graph_segment <- function(segments_data){
  plot <- ggplot() +
    geom_line(data=nb.profiles,
              mapping=aes(x=position, y=logratio), color="#636363")+
    geom_segment(data=segments_data,
                 mapping=aes(x=nb.profiles[['position']][start], y=mean,
                             xend=nb.profiles[['position']][end], yend=mean),
                 color="#fd8d3c", size=0.8) +
    geom_vline(data=segments_data[start > 1],
               mapping=aes(xintercept=nb.profiles[['position']][start]),
               color="#56B4E9", size=0.9) +
    facet_grid(segments ~ .) +
    theme_bw()
  return(plot)
}

graph_segment(models.binsegRcpp.coef)

models.changepoint <- cpt.meanvar(data=nb.profiles[['logratio']], penalty="None", pen.value=0, method="BinSeg",
                Q=4, test.stat="Normal", class=TRUE, param.estimates=TRUE)


models.changepoint.coef <- data.table(segments=integer(0), start=integer(0), end=integer(0), mean=double(0))
for (i in 1:4){
  start <- sort(na.omit(c(1, cpts.full(models.changepoint)[i, ]+1)))
  end <- sort(na.omit(c(cpts.full(models.changepoint)[i, ], length(data.set(models.changepoint)))))
  for(j in 1:(i+1)){
    tmp.mean <- mean(nb.profiles[['logratio']][start[j]:end[j]])
    models.changepoint.coef <- rbind(models.changepoint.coef, list(i+1 , start[j], end[j], tmp.mean))
  }
}

graph_segment(models.changepoint.coef)

# ggplot() +
#   geom_line(data=nb.profiles,
#             mapping=aes(x=position, y=logratio), color="#636363")+
#   geom_segment(data=segments_data,
#                mapping=aes(x=nb.profiles[['position']][start], y=mean,
#                            xend=nb.profiles[['position']][end], yend=mean),
#                color="#fd8d3c", size=0.8) +
#   geom_vline(data=segments_data[start > 1],
#              mapping=aes(xintercept=nb.profiles[['position']][start]),
#              color="#56B4E9", size=0.9) +
#   facet_grid(segments ~ .) +
#   theme_bw()
# return(plot)


binseg.cpts = data.frame(changepoint=c(1:5), index=sort(models.binsegRcpp[['end']]))
changepoint.cpts = data.frame(changepoint=c(1:4), index=cpts((models.changepoint)))

comparative_graph <- function(data1, data2, x1, y1, x2, y2, name, xlabel, ylabel){
  plot <- ggplot() +
    geom_line(data=data1,
              mapping=aes(x=data1[[x1]], y=data1[[y1]], color="binsegRcpp"), size=2) +
    geom_text(data=data1,
              aes(x=data1[[x1]]-0.2, y=data1[[y1]]+0.5,label=index),vjust=0) +
    geom_line(data=data2,
              mapping=aes(x=data2[[x2]], y=data2[[y2]], color="changepoint"), size=1) +
    geom_text(data=data2,
              aes(x=data2[[x2]]-0.2, y=data2[[y2]]+0.5,label=index),vjust=0) +
    ggtitle(label=name) +
    labs(x=xlabel, y =ylabel, colour="Package")+
    theme_bw() +
    theme(plot.title=element_text(face="bold",margin=margin(t=10,b=15),size=16),
          axis.title.x=element_text(margin=margin(t=10,b=10), size=14),
          axis.title.y=element_text(margin=margin(r=10,l=15), size=14))
    return(plot)
}

g <- comparative_graph(binseg.cpts, changepoint.cpts, 'changepoint', 'index', 'changepoint', 'index', "Predicted changepeoints per package", "Changepoint number", "Changepoint index")

changepoint.loss = data.frame(changepoint=2:5, loss=pen.value.full(models.changepoint))

ggplot() +
  geom_line(data=models.binsegRcpp,
            mapping=aes(x=segments, y = loss, color="binsegRcpp"), size=1.5) +
  geom_text(data=models.binsegRcpp,
            aes(x=segments+0.2,y=loss+0.3,label=round(loss, 2)),vjust=0) +
  geom_line(data=changepoint.loss,
            mapping=aes(x=changepoint, y = loss, color="changepoint"), size=1.5) +
  geom_text(data=changepoint.loss,
            aes(x=changepoint+0.2,y=loss+0.3,label=round(loss, 2)),vjust=0) +
  ggtitle("Loss value per segment") +
  labs(x="Segment number", y ="Loss value", colour="Package")+
  theme_bw() +
  theme(plot.title=element_text(face="bold",margin=margin(t=10,b=15),size=16),
        axis.title.x=element_text(margin=margin(t=10,b=10), size=14),
        axis.title.y=element_text(margin=margin(r=10,l=15), size=14))


negbin <- dataGenerator(n=100, changepoints=c(1), parameters=c(0.8), type="negbin", size=10)

plot(negbin)

v <- var(negbin)
mu <- mean(negbin)

p <- (mu)/v
r <- mu**2/(v-mu)

install_github("diego-urgell/binsegRcpp")
library("binsegRcpp")

# Computing the model with the modified version of binsegRcpp
models.binsegRcpp.meanvar <- binseg_normal(nb.profiles[['logratio']], 5)
models.binsegRcpp.meanvar.coef <- coef(models.binsegRcpp)

# Computing the model with the changepeoint package using change in mean and variance
models.changepoint.meanvar <- cpt.meanvar(data=nb.profiles[['logratio']], penalty="None", pen.value=0, method="BinSeg",
                                          Q=4, test.stat="Normal", class=TRUE, param.estimates=TRUE)

cpts <- ggplot() +
  geom_line(data=data.table(s = sort(models.binsegRcpp.meanvar$end)), mapping=aes(x=1:5, y=s))

cpts


print(p)
print(r)

cpts <- ggplot() +
  geom_line(data=data.table(s = sort(models.binsegRcpp.meanvar$end)),
            mapping=aes(x=1:5, y=s, color="binsegRcpp")) +
  geom_text(data=data.table(s = sort(models.binsegRcpp.meanvar$end)),
                            aes(x=(1:5)+0.2,y=s,label=round(s, 2)),vjust=0) +
  geom_line(data=data.table(s = unlist(cpts(models.changepoint.meanvar))),
            mapping=aes(x=1:4, y=s, color="changepoint")) +
  geom_text(data=data.table(s = unlist(cpts(models.changepoint.meanvar))),
            aes(x=(1:4)+0.2,y=s,label=round(s, 2)),vjust=0) +
  ggtitle("Changepoints per package") +
  labs(x="Segment number", y ="Loss value", colour="Package")+
  theme_bw() +
  theme(plot.title=element_text(face="bold",margin=margin(t=10,b=15),size=16),
        axis.title.x=element_text(margin=margin(t=10,b=10), size=14),
        axis.title.y=element_text(margin=margin(r=10,l=15), size=14))

cpts



#
#
#
# for(i in 1:4){
#   tmp <- cpt.mean(data=nb.profiles[['logratio']], penalty="None", pen.value=0, method="BinSeg",
#                   Q=i, test.stat="Normal", class=TRUE, param.estimates=TRUE)
#   start <- unlist(c(1, cpts(tmp)[1:i]))
#   end <- c(unlist(cpts(tmp)), length(data.set(tmp)))
#   mean <- unlist(param.est(tmp))
#   tmp_table <- data.table(rep(i+1, i+1), start, end, mean)
#   names(tmp_table) <- name
#   models.changepoint <- rbind(models.changepoint, tmp_table)
# }
#
#
#
# #data <- c(rnorm(100, 15), rnorm(100, 0))
# #model <- cpt.mean(data=data, penalty="None", pen.value=0, method="BinSeg", Q=1, test.stat="Normal", class=TRUE, param.estimates=TRUE)


