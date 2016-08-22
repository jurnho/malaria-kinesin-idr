disprot = read.csv("disprot.csv", stringsAsFactors=FALSE)
protein_ids = unique(disprot$protein_id)

prediction_files = Sys.glob("DP*csv")
# reduce to big matrix of predictions.

# list of dataframes
prediction_dfs = lapply(prediction_files, function(prediction_file) {
    read.csv(prediction_file, stringsAsFactors=FALSE)
})
# https://stat.ethz.ch/pipermail/r-help/2010-September/252046.html
predictions = do.call("rbind", prediction_dfs)

## output results_summary
# protein_id, predictor, tp, fn, fp, n
results_summary = data.frame(protein_id=character(), uniprot_id=character(), name=character(),
predictor=character(), TP=numeric(), FN=numeric(), FP=numeric(), TN=numeric(), sensitivity=numeric(),
specificity=numeric(), accuracy=numeric(), pdb_id=character(), disprot_disorder_count=character(),
disprot_order_count=character()
)

# threshold is not always 0.5, it depends on the predictor
is_disorder_predicted = function(predictor, raw_scores) {
    threshold = 0.5
    if (predictor == 'DISEMBL_COILS') {
        threshold = 0.43 * 1.1
    }
    if (predictor == 'DISEMBL_REM465') {
        threshold = 0.5 * 1.1
    }
    if (predictor == 'DISEMBL_HOTLOOPS') {
        threshold = 0.086 * 1.2
    }
    is_disordered = function() {
        if (raw_score >= threshold) {
            return (TRUE)
        }
        return (FALSE)
    }

    ordered = raw_scores >= threshold
    return (ordered)
}

for (protein_id in protein_ids) {
    message(protein_id)
    uniprot_id = disprot[disprot$protein_id == protein_id,][1,"uniprot_id"]
    name = disprot[disprot$protein_id == protein_id,][1,"name"]
    pdb_id = disprot[disprot$protein_id == protein_id,][1,"pdb.id"]
    disprot_disorder_count = sum(disprot[disprot$protein_id==protein_id,c('disordered')] == 'Y')
    disprot_order_count = sum(disprot[disprot$protein_id==protein_id,c('disordered')] == 'N')

    disprot_protein_info = disprot[disprot$protein_id == protein_id,]
#    message(length(disprot_protein_info$disordered))
    # predictions for a single protein with all predictors
    protein_prediction_results = predictions[predictions$protein_id == protein_id,]
    predictors = unique(protein_prediction_results$predictor)
    for (predictor in predictors) {
        message(predictor)
        # predictions for a a single protein and single predictor
        predictor_results = protein_prediction_results[protein_prediction_results$predictor == predictor,]
#        message(length(predictor_results$score))
        if (length(disprot_protein_info$disordered) != length(predictor_results$score)) {
          message("wrong count")
          next
        }
        disorder_predicted = is_disorder_predicted(predictor, predictor_results$score)
        true_positive = sum(disorder_predicted & (disprot_protein_info$disordered =='Y'))
        # false negative
        false_negative = sum((!disorder_predicted) & (disprot_protein_info$disordered =='Y'))
        # false positive
        false_positive = sum(disorder_predicted & (disprot_protein_info$disordered =='N'))
        # true negative
        true_negative = sum((!disorder_predicted) & (disprot_protein_info$disordered =='N'))

        sensitivity = true_positive / (true_positive + false_negative)
        specificity = true_negative / (true_negative + false_positive)
        accuracy = (true_positive + true_negative) / (true_positive + false_positive + false_negative + true_negative)
        #    message(paste(protein_id, predictor, true_positive, false_negative, false_positive, true_negative, sep=","))
        #    v = c(protein_id, predictor, true_positive, false_negative, false_positive, true_negative)
        #    names(v) = c('protein_id','predictor','TP','FN','FP','TN')
        results_summary = rbind(results_summary, data.frame(
        protein_id = protein_id, uniprot_id=uniprot_id, name=name, predictor = predictor, TP = true_positive, FN = false_negative,
        FP = false_positive, TN = true_negative, sensitivity = sensitivity, specificity = specificity, accuracy = accuracy,
        pdb_id = pdb_id, disprot_disorder_count = disprot_disorder_count, disprot_order_count = disprot_order_count
        ))
    }
}

write.csv(results_summary, file = "predictor_disprot_results_summary.csv")

# results_summary contains accuracy by predictor and protein

predictors = unique(results_summary$predictor)
predictor_scores = data.frame(predictor=character(), test_set=character(),
sensitivity=numeric(), specificity=numeric(), accuracy=numeric()
)

predictor_mean_score = function(protein_scores, test_set_name, predictor) {
    p_sensitivity = protein_scores$sensitivity
    p_specificity = protein_scores$specificity
    p_accuracy = protein_scores$accuracy

    mean_sensitivity = mean(p_sensitivity[!is.nan(p_sensitivity)])
    mean_specificity = mean(p_specificity[!is.nan(p_specificity)])
    mean_accuracy = mean(p_accuracy[!is.nan(p_accuracy)])
    p_result = data.frame(predictor=predictor,
        test_set=test_set_name,
        mean_sensitivity=mean_sensitivity,
        mean_specificity=mean_specificity,
        mean_accuracy=mean_accuracy
    )
    return (p_result)
}
for (predictor in predictors) {
    p = results_summary[results_summary$predictor==predictor,]
    scores = predictor_mean_score(p, 'Disprot all', predictor)
    predictor_scores = rbind(predictor_scores, scores)
    p_pdb = p[p$pdb_id!="",]
    scores = predictor_mean_score(p_pdb, 'Disprot with PDB only', predictor)
    predictor_scores = rbind(predictor_scores, scores)
}

predictor_scores

write.csv(predictor_scores, file = "predictor_scores.csv")

# correlations from results
#vector
correlations = c()
for (p1 in predictors) {
    s1 = predictions[predictions$predictor==p1,]
    for (p2 in predictors) {
        s2 = predictions[predictions$predictor==p2,]
        j3 = merge(s1,s2,by=c('protein_id','position'))
        correlation = cor(j3$score.x, j3$score.y)
        message(p1)
        message(p2)
        message(correlation)
        correlations = c(correlations, correlation)
    }
}
lpred = length(predictors)
m = matrix(nrow=lpred, ncol=lpred, data=correlations)
rownames(m) = predictors;
colnames(m) = predictors;
write.csv(m, file = "predictor_correlations.csv")

disprot_predictor_scores = disprot
for (p in predictors) {
    scores = predictions[predictions$predictor==p,c('protein_id', 'position','score')]
    disprot_predictor_scores = merge(disprot_predictor_scores,scores,by=c('protein_id','position'))
    # update column names with name of predictor
    new_column_names = c(head(colnames(disprot_predictor_scores),-1), p)
    colnames(disprot_predictor_scores) = new_column_names
}

write.csv(disprot_predictor_scores,file='disprot_predictor_scores.csv')

# only keep proteins with a PDB, these are more simliar to kinesins
pdb_only=disprot_predictor_scores[disprot_predictor_scores$pdb.id!='',]
# disordered residue count
sum(pdb_only$disordered=='Y')
# ordered residue count
sum(pdb_only$disordered=='N')

protein_ids = unique(pdb_only$protein_id)
# randomly select 70% of proteins for training, and 30% for testing
buckets = sample(2, length(protein_ids), replace=TRUE, prob=c(0.70, 0.30))
protein_buckets = cbind(protein_ids, buckets)
colnames(protein_buckets) = c('protein_id','bucket')

# run knn on the disprot dataset
knn_disprot = function(pdb_only, protein_buckets, k=3) {
  knn_data_set = merge(pdb_only,protein_buckets,by=c('protein_id'))
  library(class)
  knn_training = knn_data_set[knn_data_set$bucket==1,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
  knn_training.labels = knn_data_set[knn_data_set$bucket==1,c('disordered')]
  knn_testing = knn_data_set[knn_data_set$bucket==2,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
  knn_testing.labels = knn_data_set[knn_data_set$bucket==2,c('disordered')]
  knn_result = knn(train=knn_training, test=knn_testing, cl=knn_training.labels, k)

  true_positive = sum((knn_result=='Y') & (knn_testing.labels =='Y'))
  # false negative
  false_negative = sum((knn_result=='N') & (knn_testing.labels =='Y'))
  # false positive
  false_positive = sum((knn_result =='Y') & (knn_testing.labels =='N'))
  # true negative
  true_negative = sum((knn_result =='N') & (knn_testing.labels =='N'))

  knn_sensitivity = true_positive / (true_positive + false_negative)
  knn_specificity = true_negative / (true_negative + false_positive)
  knn_accuracy = (true_positive + true_negative) / (true_positive + false_positive + false_negative + true_negative)
  result = data.frame(k=k,
    knn_sensitivity=round(knn_sensitivity, 3),
    knn_specificity=round(knn_specificity, 3),
    knn_accuracy=round(knn_accuracy,3)
  )
  return(result)
}
knn_disprot(pdb_only, protein_buckets, 3)
knn_disprot(pdb_only, protein_buckets, 5)
knn_disprot(pdb_only, protein_buckets, 7)
knn_disprot(pdb_only, protein_buckets, 9)

y_score = function(disordered_YN) {
    if (disordered_YN == 'Y') {
        return (1)
    }
    return (0)
}

to_disorder_YN = function(lm_score) {
    if (lm_score > 0.5) {
        return ('Y')
    }
    return ('N')
}

# linear regression
lm_disprot = function(pdb_only, protein_buckets) {
  lm_y = sapply(pdb_only$disordered, y_score)
  lm_y = as.vector(lm_y)
  # add a new y_score column as 1 or disorder, 0 for ordered.
  disprot_lm = cbind(pdb_only, lm_y)
  lm_data_set = merge(disprot_lm,protein_buckets,by=c('protein_id'))
  lm_data = lm_data_set[lm_data_set$bucket==1,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short','lm_y')]
  # use all of the above
  lm_fit = lm(lm_y ~ ., data = lm_data)
  message(coef(lm_fit))
  lm_test_data = lm_data_set[lm_data_set$bucket==2,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
  lm_test_labels = lm_data_set[lm_data_set$bucket==2,c('disordered')]
  lm_predictions = predict(lm_fit, lm_test_data)
  lm_predictions_YN = sapply(lm_predictions,to_disorder_YN)
  lm_predictions_YN = as.vector(lm_predictions_YN)
  lm_results = lm_data_set[lm_data_set$bucket==2,]
  lm_results = cbind(lm_results, lm_predictions, lm_predictions_YN)
  lm_results = lm_results[order(lm_results["protein_id"], lm_results["position"]),]
  write.csv(lm_results, file = "lm_results.csv")
  true_positive = sum((lm_predictions_YN=='Y') & (lm_test_labels =='Y'))
  # false negative
  false_negative = sum((lm_predictions_YN=='N') & (lm_test_labels =='Y'))
  # false positive
  false_positive = sum((lm_predictions_YN =='Y') & (lm_test_labels =='N'))
  # true negative
  true_negative = sum((lm_predictions_YN =='N') & (lm_test_labels =='N'))
  lm_sensitivity = true_positive / (true_positive + false_negative)
  lm_specificity = true_negative / (true_negative + false_positive)
  lm_accuracy = (true_positive + true_negative) / (true_positive + false_positive + false_negative + true_negative)
  result = data.frame(
    lm_sensitivity = round(lm_sensitivity, 3),
    lm_specificity = round(lm_specificity, 3),
    lm_accuracy = round(lm_accuracy, 3)
  )
  return (result)
}
lm_disprot(pdb_only, protein_buckets)

# results is a vector of Y or N
# remove lone residue predictions of disorder and order.
smooth_regions = function(predictions) {
  # loop through each residue.
  lone_disorder_removed = c()
  for (i in seq(1, length(predictions))) {
    if (i==1 | i==length(predictions)) {
      lone_disorder_removed = c(lone_disorder_removed, predictions[i])
    }
    # else look at +1 either side
    else {
      window = c(predictions[i-1], predictions[i], predictions[i+1])
      # if it is inbetween disordered regions, set it as disordered
      if (identical(window, c('N','Y','N'))) {
        lone_disorder_removed = c(lone_disorder_removed,'N')
      }
      else {
        lone_disorder_removed = c(lone_disorder_removed, predictions[i])
      }
    }
  }
  # again but this time for ordered residues (N)
  results = c()
  for (i in seq(1, length(lone_disorder_removed))) {
    if (i==1 | i==length(lone_disorder_removed)) {
      results = c(results, lone_disorder_removed[i])
    }
    # else look at +1 either side
    else {
      window = c(lone_disorder_removed[i-1], lone_disorder_removed[i], lone_disorder_removed[i+1])
      # if it is inbetween disordered regions, set it as disordered
      if (identical(window, c('Y','N','Y'))) {
        results = c(results,'N')
      }
      else {
        results = c(results, lone_disorder_removed[i])
      }
    }
  }
  return (results)
}


# create a SVM and evaluate the results.
svm_disprot = function(pdb_only, protein_buckets, yes_weight=4, num_samples=5000, svm_kernel="radial",
training_subset=NULL, filter_lone_regions=TRUE) {
  library("e1071")
  svm_data_set= merge(pdb_only,protein_buckets,by=c('protein_id'))
  svm_train.all = svm_data_set[svm_data_set$bucket==1,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short','disordered')]
  # random sample from the training subset
  if (is.null(training_subset)) {
    training_subset = sample(nrow(svm_train.all), num_samples)
  }
  svm_train.subset = svm_train.all[training_subset,]
  svm_train.data = svm_train.subset[,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
  svm_train.labels = as.factor(svm_train.subset[,c('disordered')])
  svm_model = svm(x=svm_train.data, y=svm_train.labels, probability = TRUE, class.weights=c('Y'=yes_weight, 'N'=1),
                kernel=svm_kernel)
  svm_test_all = svm_data_set[svm_data_set$bucket==2,c('protein_id','position','disordered','DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
  svm_test_all = svm_test_all[order(svm_test_all["protein_id"], svm_test_all["position"]),]
  svm_test_data = svm_test_all[,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
  svm_test_labels = svm_test_all[,c('disordered')]
  # probability=TRUE doesn't seem to work here.
  svm_prediction = predict(svm_model, svm_test_data, decision.values = TRUE)
  svm_prediction_YN = as.vector(svm_prediction)
  svm_prediction_decision_values = attr(svm_prediction, "decision.values")[,1]
  # sometimes the class label is reversed
  if (colnames(attr(svm_prediction, "decision.values"))[1] == 'N/Y') {
      svm_prediction_decision_values = svm_prediction_decision_values * -1
  }
  svm_results = cbind(svm_test_all, svm_prediction_YN)
  if (filter_lone_regions) {
    svm_prediction_YN = smooth_regions(svm_prediction_YN);
  }
  # filter svm_prediction_YN
  true_positive = sum((svm_prediction_YN=='Y') & (svm_test_labels =='Y'))
  false_negative = sum((svm_prediction_YN=='N') & (svm_test_labels =='Y'))
  false_positive = sum((svm_prediction_YN =='Y') & (svm_test_labels =='N'))
  true_negative = sum((svm_prediction_YN =='N') & (svm_test_labels =='N'))

  svm_sensitivity = true_positive / (true_positive + false_negative)
  svm_specificity = true_negative / (true_negative + false_positive)
  svm_accuracy = (true_positive + true_negative) / (true_positive + false_positive + false_negative + true_negative)

  result = c('svm_sensitivity'=round(svm_sensitivity, 3), 'svm_specificity'=round(svm_specificity,3),
        'svm_accuracy'=round(svm_accuracy,3))
  return (list(result=result, svm_model=svm_model,
  svm_prediction_decision_values=svm_prediction_decision_values,
  svm_prediction_YN = svm_prediction_YN,
  training_subset = training_subset))
}

svm_disprot_result = svm_disprot(pdb_only, protein_buckets, yes_weight = 4.5)
svm_result=svm_disprot_result$result
svm_model=svm_disprot_result$svm_model
svm_prediction_decision_values=svm_disprot_result$svm_prediction_decision_values
svm_prediction_YN=svm_disprot_result$svm_prediction_YN

# plasmodb kinesin files are found in plasmodb/csv/
plasmodb_kinesins = Sys.glob('plasmodb/csv/*.csv')
# many dataframes
plasmodb_ref_dfs = lapply(plasmodb_kinesins, function(plasmo) {
    read.csv(plasmo, stringsAsFactors=FALSE)
})
# combine them into one
plasmodb = do.call("rbind", plasmodb_ref_dfs)

# load up individual predictions

plasmo_prediction_files = Sys.glob('plasmodb/predictions/*/*.csv')
plasmo_prediction_dfs = lapply(plasmo_prediction_files, function(plasmo_prediction_file) {
    read.csv(plasmo_prediction_file, stringsAsFactors=FALSE)
})
plasmo_prediction_df = do.call("rbind", plasmo_prediction_dfs)

plasmodb_predictions = plasmodb
for (p in predictors) {
    # extract scores for the predictor
    scores = plasmo_prediction_df[plasmo_prediction_df$predictor==p,c('protein_id', 'position','score')]
    # join
    plasmodb_predictions = merge(plasmodb_predictions,scores,by=c('protein_id','position'))
    # update col name
    new_column_names = c(head(colnames(plasmodb_predictions),-1), p)
    colnames(plasmodb_predictions) = new_column_names
}

# svm input
plasmo_svm.test_data = plasmodb_predictions[,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
# predict
plasmo_svm.prediction = predict(svm_model, plasmo_svm.test_data, decision.values = TRUE)
# results
plasmo_svm.disorder = as.vector(plasmo_svm.prediction)

plasmo_svm_decision_values = attr(plasmo_svm.prediction, "decision.values")[,'N/Y'] * -1
plasmo_results = cbind(plasmodb_predictions, plasmo_svm.disorder)
plasmo_results = cbind(plasmo_results, plasmo_svm_decision_values)
plasmo_results = plasmo_results[order(plasmo_results["protein_id"], plasmo_results["position"]),]
write.csv(plasmo_results, file = "plasmo_svm_prediction.csv")

# other info
## percent of predicted disorder in test set
sum(svm_prediction_decision_values > 0)
sum(svm_prediction_decision_values < 0)
test_set_predicted_percent_disorder = sum(svm_prediction_decision_values > 0) / length(svm_prediction_decision_values)
round(test_set_predicted_percent_disorder,3)
## percent of predicted disorder in disprot
sum(plasmo_svm_decision_values > 0)
sum(plasmo_svm_decision_values < 0)
plasmo_predicted_percent_disorder = sum(plasmo_svm_decision_values > 0) / length(plasmo_svm_decision_values)
round(plasmo_predicted_percent_disorder,3)

# plots
breaks = seq(-2.4, 2.4, by=0.10 )
h = hist(svm_prediction_decision_values, breaks=breaks)
h$counts = h$counts / length(svm_prediction_decision_values)

h2 =hist(plasmo_svm_decision_values, breaks=breaks)
h2$counts = h2$counts / length(plasmo_svm_decision_values)

a = rbind(h$counts, h2$counts)
rownames(a) = c('test set','plasmodb kinesins')
barplot(a, beside=TRUE, xlab = "decision value (distance)", legend=rownames(a),
main="SVM relative frequency distribution of decision values", names=h$mids)
png("svm_relative_frequency.png", width = 800, height = 600)
barplot(a, beside=TRUE, xlab = "decision value (distance)", legend=rownames(a),
main="SVM relative frequency distribution of decision values", names=h$mids)
dev.off()

png(filename="predictor_score_distributions.png")
library(RColorBrewer)
colors = brewer.pal(6, "Dark2")
#colors = c("red", "blue", "green", "purple", "brown", "orange")
line_types = c(1,2,3,4,5,6)
plot(density(pdb_only$DISOPRED, adjust=0.2), col=colors[1], lty=line_types[1], lwd=2,
main = "Score distribution of individual predictors", xlab='predictor score')
lines(density(pdb_only$DISEMBL_COILS ), col=colors[2], lty=line_types[2], lwd=2)
lines(density(pdb_only$DISEMBL_REM465 ), col=colors[3], lty=line_types[3], lwd=2)
lines(density(pdb_only$DISEMBL_HOTLOOPS ), col=colors[4], lty=line_types[4], lwd=2)
lines(density(pdb_only$iupred_long), col=colors[5], lty=line_types[5], lwd=2)
lines(density(pdb_only$iupred_short), col=colors[6], lty=line_types[6], lwd=2)
legend(x='topright', legend = c('DISOPRED','DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','iupred_long','iupred_short'), lty=line_types, col=colors, lwd=2)
dev.off()

# return svm roc data for plotting
svm_roc = function(kernel='radial', training_subset = NULL) {
  # we can tweak the weight and see the ROC changes.
  yes_weights=c(1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,9,10,15,20)
  sensitivities=c()
  specificities=c()
  for (yes_weight in yes_weights) {
    s = svm_disprot(pdb_only, protein_buckets, yes_weight=yes_weight, num_samples=5000, svm_kernel=kernel,
                    training_subset=training_subset)
    sensitivities = c(sensitivities, s$result['svm_sensitivity'])
    specificities = c(specificities, s$result['svm_specificity'])
  }
  return (list(sensitivities=as.vector(sensitivities), specificities=as.vector(specificities), yes_weights=yes_weights))
}

predictor_test_data_scores = function() {
  test_data_set = merge(pdb_only,protein_buckets,by=c('protein_id'))
  test_data_set = test_data_set[test_data_set$bucket==2,c('protein_id','position','disordered','DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
  sensitivities=c()
  specificities=c()
  for (predictor in predictors) {
    predictions = is_disorder_predicted(predictor, test_data_set[,c(predictor)])
    actuals = (test_data_set$disordered == 'Y')
    true_positive = sum(predictions & actuals)
    false_positive = sum(predictions & !actuals)
    true_negative = sum(!predictions & !actuals)
    false_negative = sum(!predictions & actuals)
    sensitivity = true_positive / (true_positive + false_negative)
    specificity = true_negative / (true_negative + false_positive)
    sensitivities = c(sensitivities, sensitivity)
    specificities = c(specificities, specificity)
  }
  return (list(sensitivities=sensitivities, specificities = specificities, predictors = as.vector(predictors)))
}

# roc for default radial SVM
roc_predictor = predictor_test_data_scores()
roc = svm_roc()

png(filename="svm_roc.png", width = 800, height = 600)
plot(1-roc$specificities, roc$sensitivities, xlim=c(0,1), ylim=c(0,1),
main="ROC curve for SVM trained with various disorder class weights", type="l", col="blue",ylab='sensitivity',
xlab='1-specificity')
points(1-roc$specificities, roc$sensitivities,pch=1, col="blue")
standalone_sensitivities = roc_predictor$sensitivities
standalone_specificities = roc_predictor$specificities
# pch=4 cross
points(1-standalone_specificities, standalone_sensitivities, pch=4, col="darkgreen")
abline(0, 1,col="red")
legend(x='bottomright', legend = c('SVM consensus predictor','individual predictors','random predictor'),
col=c('blue','darkgreen','red'), lty=c(1,NA,1), lwd=2, pch=c(1,4,NA))
dev.off()

#x11()

# run it once and remember the training subset
ts = svm_disprot(pdb_only, protein_buckets)$training_subset

radial_roc = svm_roc(kernel='radial', training_subset =ts)
polynomial_roc = svm_roc(kernel='polynomial', training_subset =ts)
linear_roc = svm_roc(kernel='linear', training_subset =ts)
sigmoid_roc = svm_roc(kernel='sigmoid', training_subset =ts)
colors = brewer.pal(4, "Dark2")
line_types = c(1,2,3,4)

png(filename="svm_roc_kernels.png", width = 800, height = 600)

plot(1-radial_roc$specificities, radial_roc$sensitivities, xlim=c(0,1), ylim=c(0,1),
main="ROC curve for SVM trained with various kernels", type="l", col=colors[1],ylab='sensitivity',
xlab='1-specificity', lty=line_types[1],lwd=2)
lines(1-polynomial_roc$specificities, polynomial_roc$sensitivities, col=colors[2],lty=line_types[2],lwd=2)
lines(1-linear_roc$specificities, linear_roc$sensitivities, col=colors[3], lty=line_types[3],lwd=2)
lines(1-sigmoid_roc$specificities, sigmoid_roc$sensitivities, col=colors[4], lty=line_types[4],lwd=2)
legend(x='bottomright', legend = c('radial','polynomial','linear','sigmoid'),
 col=colors, lty = line_types, lwd=2)
dev.off()

# check difference between filter and no filter
no_filter = svm_disprot(pdb_only, protein_buckets, yes_weight=4.5, filter_lone_regions=FALSE)
svm_disprot(pdb_only, protein_buckets, yes_weight=4.5, training_subset=no_filter$training_subset)$result

x11()