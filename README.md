# Peril-incendie
## Résumé
Le péril incendie correspond à un péril majeur en assurance dommages. Il est une des raisons qui ont poussé au développement de l'assurance de dommage. Dans ce rapport, on utilise plusieurs familles paramétriques pour modéliser la sévérité des données canadiennes de sinistre incendie. Les lois de probabilité considérées sont la lognormale-Pareto généralisée, Weibull-Pareto généralisée, Coxienne-2-Pareto généralisée et la Bêta généralisée de type 2-Pareto généralisée. On utilise le maximum de vraisemblance pour estimer les paramètres des lois et l’on utilise des tests d'adéquations graphiques en plus de tests quantitatifs pour juger de la performance des lois ajustées aux données. On calcule les mesures $VaR$ et $TVaR$ pour avoir de l'information sur la queue de la distribution. Pour le processus d'avènement des incendies, on utilise un processus de Poisson homogène pour la modélisation. Après avoir trouvé une loi pour le processus d'avènement des sinistres et une loi pour le montant des sinistres, on veut créer un fonds fictif pour analyser les entrées et sorties de fonds. À l'aide de la méthode de simulation Monte-Carlo, on peut approximer la probabilité de ruine ce qui permettra d'avoir une idée des chances que le fonds devienne insolvable.

## Description des documents 
Une brève description des documents qui compose le dépot
- Code : toutes les fichiers des code neccessaire pour le projet
    - script.R : Le code principale du projet
- Article : article lu et utilisé pour le projet
- Data : les bases de données utiliser pour le projet
- Image : les graphiques produit pour le projet
