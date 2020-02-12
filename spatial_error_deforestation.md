# Performance evaluation for forecasting modeling with spatiotemporal structures in data
Yujun Zhou and Kathy Baylis

## Research question

What are the preferred cross-validation (CV) and out-of-sample (OOS) performance evaluation methods for given spatiotemporal correlations in panel data?

How to adjust for spatiotemporal correlations to improve the out-of-sample performance of the models?

How to tackle spatial and temporal heterogeneity using clustering?

## Motivation

Geospatial data have become widely available in recent years and have made possible more spatial-temporal forecasting applications, including deforestation prediction and air pollution monitoring. Choosing the right performance evaluation matters for generating accurate and trustworthy predictions. Two primary evaluation methods are cross-validation and Out-of-sample evaluation. Cross-validation (CV) is a resampling-based technique for the estimation of a model’s predictive performance. K-fold CV, for example, divide the data into K subsets of approximately the same size and then having each subset used successively as the test set. Out-of-sample (OOS) evaluation is to estimate the performance of the model in “unseen” data. The standard way of doing this is called “hold-out validation,” namely splitting the data into training for learning and untouched testing set hold out for evaluation. With spatial-temporal dependences between observations in both training and testing data, model performance evaluated using CV and OOS can be problematic. In both cases, the procedures random resampling of the training data or random splitting into the testing data are no longer random. In other words, the observations in the testing set are no longer independent from the training set with the spatiotemporal correlations (Oliveira et al., 2018).

As is described in previous works (particularly Robert et al., 2017), the consequence of ignoring these correlations include the following:

1. Unreliable results: compared to the result tested on a real independent set, CV and OOS results would appear to be better than they are, making predictions less reliable.

2. Overfitting: more complicated models are preferred in the context of the correlations. Better performance on the training are closely associated with better performance on the testing set, more so because the independence assumption is violated.

3. Misinterpretation of variables: The correlations can be incorporated into the model process through some of the covariates that also have spatial-temporal variation and appear to be more predictive than they actually are, e.g. rainfall.

However, the severity of consequences is not apparent for different spatial and temporal structures. In this study, we address this issue using Monte Carlo simulation to specify different correlation structures and varying amounts of correlations between observations. This allows us to see how significant bias can be across different cases. Bearing the same correlation structures in mind, we assess the results on real-world datasets in household surveys and deforestation. The performance is compared to a “truly” independent dataset created using stratified random sample of the original test set with the proposed correlations.
Previous studies have come up with several ways of addressing spatial and temporal autocorrelations with similar ideas. Brenning (2012) uses k-means clustering to partition to reduce the influence of spatial autocorrelation. Roberts et al. (2017) invent “block cross-validation” by introducing spatial blocks on contiguous geographic space to force the model to be tested on more distant data, with similar ideas for temporal blocks of grouping data that fall in the nearby time interval. Meyer et al. (2018) create “Leave-Location-and-Time-Out (LLTO) CV” to address spatiotemporal correlations. Schratz et al. (2019) focus on the hyperparameter tuning aspect concerning the underlying spatial correlations in the data. We compare the effectiveness of these methods in both the simulated dataset and the real-world panel data.

One closely related but under-addressed issue is the heterogeneity in the panel data. Such heterogeneity can come from seasonality (e.g., air pollution in winter), natural disasters (severe drought in one year), geography (flood in low altitude regions). All of these cases would make it hard to forecast on the test dataset. Making use of the spatial and temporal correlations within the blocks can help improve model performance. Specifically, we can do so by taking groups of observations that are correlated as 'clusters' in the train set with strong spatial correlation, or years together for the same observation if there are strong dynamic effects.

## Contribution

Explore different spatial-temporal error structures in panel dataset how they affect the model performance, model selection differently

Showcase the bias and apply the adjustment for spatial-temporal correlation in both environmental and development settings in real-world data

Adjustment for spatiotemporal heterogeneity in the data

## Methodology

###	Monte Carlo simulation

Consider a Data Generating Process (DGP) in a panel setting with spatial-temporal correlations:

$$Y_{it}\ =  \beta\ X_{it}\ +\ u_i\ +\ v_t\ +\ \gamma_{it}{+\ \varepsilon}_{it} $$

where u_i follows a standard spatial lag process, v_t follows a standard temporal lag, $\gamma_{it}$ follow a spatial diffusion process (correlated with $\gamma_{it-1}$) and a true random error ${\varepsilon}_{it}$.

4.1.1. Define different spatial-temporal error structures (allowing for spatial diffusion patterns, in the case of air pollution among others)

Deforestation patterns.

	1.	Fish-bone.

Deforestation due to build or paved roads that generate side roads used by loggers or miners to reach their quarries, consequently this side roads spawn smaller ones creating the deforestation pattern of a “fish bone” (Figure 1).  95% of deforestation in Brazil Amazon has taken place within five and a half kilometers of a road or one kilometer of a navigable river. (F. Seymour, J. Busch. Why Forest, Why now, 2016)

	2.	Radial.

Radial (or pie) patterns of land use is typical of new settlements or resettlement schemes and large agricultural development. For example, San Javier in the Santa Cruz Department of Bolivia, where tropical forest was cleared in the mid- 80’s to resettle people from the Altiplano in the context of the Tierra Baja project. At the center of each unit is a small community including a church, bar/café, school and soccer field (essentials of life in rural Bolivia)


Figure 2, Source: Earth Observatory, Deforestation in Bolivia


	3.	Dendritic

Related to deforestation caused by miners and loggers that pushed into undeveloped parts of the forest, building roads that generally followed the natural contours of the land and formed curving, dendritic shapes.

Figure 3, Source: Earth Observatory, Making Sense of Amazon Deforestation Patterns

	4.	Rectangular

This pattern of deforestation tends to leave large rectangular clearings that reflect careful surveying by large-scale cattle-ranching operations. In Paraguay deforestation has been widespread in the last years. From 1987 to 2012 Paraguay has lost near 44,000 square kilometers (17,000 square miles), mainly because of the expansion of cattle farms in the western part of the country.

4.1.2. Demonstration of consequences for ignoring the spatial-temporal correlation

4.1.3. Adjustment for spatial-temporal error structures (k-means clustering, blocking)

4.1.4. Clustering methods for dealing with heterogeneity


###	Empirical data application
Evaluate on two types of widely used panel datasets in development and environmental studies: household surveys and deforestation.

Assess spatial-temporal correlations in the data

4.2.2    CV and OOS performance comparison (on a correlated set and a relatively independent set)

4.2.3.   model implications (model complexity, variable selection, the spatial distribution of prediction error)

4.2.4.   adjustment for possible spatial-temporal correlations (k-means clustering, blocking)

4.2.5.   clustering methods for dealing with heterogeneity

## Data

Simulated data with known data generated process on the spatial-temporal correlation

LSMS (living standards monitoring survey) and Global Forest Change data (2000-2017)


## References

Brenning, Alexander. "Spatial cross-validation and bootstrap for the assessment of prediction rules in remote sensing: The R package sperrorest." In 2012 IEEE international geoscience and remote sensing symposium, pp. 5372-5375. IEEE, 2012.

Meyer, Hanna, Christoph Reudenbach, Tomislav Hengl, Marwan Katurji, and Thomas Nauss. "Improving performance of spatio-temporal machine learning models using forward feature selection and target-oriented validation." Environmental Modelling & Software 101 (2018): 1-9.

Roberts, David R., Volker Bahn, Simone Ciuti, Mark S. Boyce, Jane Elith, Gurutzeta Guillera‐Arroita, Severin Hauenstein, et al. "Cross‐validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure." Ecography 40, no. 8 (2017): 913-929.

Oliveira, Mariana, Luís Torgo, and Vítor Santos Costa. "Evaluation procedures for forecasting with spatio-temporal data." In Joint European Conference on Machine Learning and Knowledge Discovery in Databases, pp. 703-718. Springer, Cham, 2018.

Schratz, Patrick, Jannes Muenchow, Eugenia Iturritxa, Jakob Richter, and Alexander Brenning. "Performance evaluation and hyperparameter tuning of statistical and machine-learning models using spatial data." arXiv preprint arXiv:1803.11266 (2018).
