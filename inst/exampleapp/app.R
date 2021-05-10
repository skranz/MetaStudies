library(MetaStudies)

# We the MetaStudy App pre-loading an example csv file
# We also add the specification test tab that shows different correlations
# between the logs of the estimates and standard errors attempting
# to gain information whether both variables are indeed independent
# in the latent distributionn absent publication bias as the Andrews & Kasy (2019)
# approach requires.

example.csv = system.file("data/IV.csv", package="MetaStudies")
MetaStudiesApp(example.csv, show.cor=TRUE)

# To call the app in its original form, you could call
#
# MetaStudiesApp(show.cor=FALSE)
