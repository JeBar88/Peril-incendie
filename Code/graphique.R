str(don$TFS_Alarm_Time)
patron <- "%Y-%m-%dT%H:%M:%S"
don$TFS_Alarm_Time <- as.Date(don$TFS_Alarm_Time, patron, tz = "UTC")
don$years <- as.factor(format(don$TFS_Alarm_Time, format = "%Y"))
perte <- don$Estimated_Dollar_Loss[!is.na(don$Estimated_Dollar_Loss) & don$Estimated_Dollar_Loss > 0]
Annee <- don$years[!is.na(don$Estimated_Dollar_Loss) & don$Estimated_Dollar_Loss > 0]
date <- don$TFS_Alarm_Time[!is.na(don$Estimated_Dollar_Loss) & don$Estimated_Dollar_Loss > 0]
status <- don$Building_Status[!is.na(don$Estimated_Dollar_Loss) & don$Estimated_Dollar_Loss > 0]


X <- data.frame(perte, date, Annee, status)



#########################################################################

Int <- c(0, 5e4, 5e5, 1e6, 1e7, 5e7, Inf)
intervalle <- function(int){
  tab <- matrix(numeric(0), ncol = 9, nrow = (length(int) - 1))
  for (i in (seq_along(int) -1)) {
   tab[i,] <- tapply(X$perte, X$Annee, function(x) length(x[x >= int[i] & x < int[i+1]]))
  }
  tab
}

apply(intervalle(Int), 1, sum)

#########################################################################


ggplot2::ggplot(X, aes(x=Annee, y=perte))+
  geom_col() + 
  labs(x="Année", y="Perte total") +
  theme_bw()

######################## Contribution #################################################


vk <- seq(0.01, 0.1, by = 0.01)
g <- sapply(vk, function(i) sum(sort(perte, dec = T)[1:(floor(i*length(perte)))])/sum(perte))
contrib <- data.frame(vk, g)
ggplot2::ggplot(contrib, aes(x=vk, y=g)) + 
  geom_col() + 
  labs(y="Contribution au total des pertes estimées", x="Pourcentage des données",
       title="") +
  theme_bw()

####################### Pertes estimées ##################################################


borne <- c(0, 6e7)
plot(X$date[X$perte >= borne[1] & X$perte < borne[2]], 
     y=X$perte[X$perte >= borne[1] & X$perte < borne[2]], "h", 
     ylab = "Perte estimée", xlab = "Date", 
     main="Pertes estimées selon la date du sinistre")
abline(h=quantile(perte, 0.999), col= "red")
grid(col = "grey")

####################### Type de batiment ##################################################
levels(X$status)[1] <- "09 - Undetermined"

borne <- c(0, 3e6)
plot(X$date[X$perte >= borne[1] & X$perte < borne[2] & X$status == "03 - Under Construction"], 
     y=X$perte[X$perte >= borne[1] & X$perte < borne[2] & X$status == "03 - Under Construction"], 
     "h", ylab = "Perte estimée", xlab = "Date", 
     main="Pertes estimées pour les batiment en constructin selon la date du sinistre")
abline(h=quantile(perte, 0.999), col= "red")
grid(col = "grey")


##############################  Autres graphiques ###########################################


ggplot2::ggplot(don, aes(x=Status_of_Fire_On_Arrival, y=Estimated_Dollar_Loss)) + 
  geom_bar(stat = "summary", fun = "mean")+ 
  theme_bw() +
  theme(axis.text.x=element_text(angle=30,hjust=1)) + 
  labs(x="État du feu à l'arrivé", y="Perte moyenne ($)",
       title = "Perte moyenne selon l'état du feu à l'arrivé") 

ggplot(don, aes(x=Smoke_Alarm_at_Fire_Origin_Alarm_Type, y=Estimated_Dollar_Loss)) + 
  geom_bar(stat = "summary", fun = "mean")+ 
  theme_bw() +
  theme(axis.text.x=element_text(angle=20,hjust=1)) +
  labs(x="Type de détecteur de fumée", y="Perte moyenne ($)",
       title = "Perte moyenne selon le type de détecteur de fumée")

ggplot(don, aes(x=don$Sprinkler_System_Presence, y=Estimated_Dollar_Loss)) + 
  geom_bar(stat = "summary", fun = "mean")+ 
  theme_bw() + 
  theme(axis.text.x=element_text(angle=20,hjust=1)) +
  labs(x="Présence de gicleur", y="Perte moyenne ($)",
       title = "Perte moyenne selon la présence de gicleur")








