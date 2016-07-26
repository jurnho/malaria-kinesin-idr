disprot = read.csv("disprot.csv", stringsAsFactors=FALSE)
protein_ids = unique(disprot$protein_id)

prediction_files = Sys.glob("DP*csv")
# reduce to big matrix of predictions.

# list of dataframs
bbb = lapply(prediction_files, function(prediction_file) {
    read.csv(prediction_file, stringsAsFactors=FALSE)
})
# https://stat.ethz.ch/pipermail/r-help/2010-September/252046.html
all_results = do.call("rbind", bbb)

## output results_summary
# protein_id, predictor, tp, fn, fp, n
results_summary = data.frame(protein_id=character(), uniprot_id=character(), name=character(),
predictor=character(), TP=numeric(), FN=numeric(), FP=numeric(), TN=numeric(), sensitivity=numeric(),
specificity=numeric(), accuracy=numeric(), pdb_id=character(), actual_disorder_count=character(),
actual_order_count=character()
)

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
    actual_disorder_count = sum(disprot[disprot$protein_id==protein_id,c('disordered')] == 'Y')
    actual_order_count = sum(disprot[disprot$protein_id==protein_id,c('disordered')] == 'N')

    disprot_actual = disprot[disprot$protein_id == protein_id,]
    prediction_results = all_results[all_results$protein_id == protein_id,]
    predictors = unique(prediction_results$predictor)
    for (predictor in predictors) {
        predictor_results = prediction_results[prediction_results$predictor == predictor,]
        #    message("predictor_results")
        #    str(predictor_results)
        #    message("actual")
        #    str(disprot_actual$disordered)
        disorder_predicted = is_disorder_predicted(predictor, predictor_results$score)
        true_positive = sum(disorder_predicted & (disprot_actual$disordered =='Y'))
        # false negative
        false_negative = sum((!disorder_predicted) & (disprot_actual$disordered =='Y'))
        # false positive
        false_positive = sum(disorder_predicted & (disprot_actual$disordered =='N'))
        # true negative
        true_negative = sum((!disorder_predicted) & (disprot_actual$disordered =='N'))

        message(predictor)
        sensitivity = true_positive / (true_positive + false_negative)
        specificity = true_negative / (true_negative + false_positive)
        accuracy = (true_positive + true_negative) / (true_positive + false_positive + false_negative + true_negative)
        #    message(paste(protein_id, predictor, true_positive, false_negative, false_positive, true_negative, sep=","))
        #    v = c(protein_id, predictor, true_positive, false_negative, false_positive, true_negative)
        #    names(v) = c('protein_id','predictor','TP','FN','FP','TN')
        results_summary = rbind(results_summary, data.frame(
        protein_id = protein_id, uniprot_id=uniprot_id, name=name, predictor = predictor, TP = true_positive, FN = false_negative,
        FP = false_positive, TN = true_negative, sensitivity = sensitivity, specificity = specificity, accuracy = accuracy,
        pdb_id = pdb_id, actual_disorder_count = actual_disorder_count, actual_order_count = actual_order_count
        ))
    }
}

write.csv(results_summary, file = "results_summary.csv")
# head(results_summary)
# mean accuracy by predictor

predictors = unique(results_summary$predictor)
predictor_summary_scores = data.frame(predictor=character(),
sensitivity=numeric(), specificity=numeric(), accuracy=numeric(),
sensitivity.pdb=numeric(), specificity.pdb=numeric(), accuracy.pdb=numeric()
)

for (predictor in predictors) {
    p = results_summary[results_summary$predictor==predictor,]
    p_sensitivity = p$sensitivity
    p_specificity = p$specificity
    p_accuracy = p$accuracy

    mean_sensitivity = mean(p_sensitivity[!is.nan(p_sensitivity)])
    mean_specificity = mean(p_specificity[!is.nan(p_specificity)])
    mean_accuracy = mean(p_accuracy[!is.nan(p_accuracy)])
    message(predictor)
    message(mean_sensitivity)
    message(mean_specificity)
    message(mean_accuracy)
    # with pdb only
    p_pdb = p[p$pdb_id!="",]
    p_pdb_sensitivity = p_pdb$sensitivity
    p_pdb_specificity = p_pdb$specificity
    p_pdb_accuracy = p_pdb$accuracy
    mean_pdb_sensitivity = mean(p_pdb_sensitivity[!is.nan(p_pdb_sensitivity)])
    mean_pdb_specificity = mean(p_pdb_specificity[!is.nan(p_pdb_specificity)])
    mean_pdb_accuracy = mean(p_pdb_accuracy[!is.nan(p_pdb_accuracy)])
    message(mean_pdb_sensitivity)
    message(mean_pdb_specificity)
    message(mean_pdb_accuracy)
    # add to results
    predictor_summary_scores = rbind(predictor_summary_scores,
    data.frame(predictor=predictor,
    mean_sensitivity=mean_sensitivity,
    mean_specificity=mean_specificity,
    mean_accuracy=mean_accuracy,
    mean_pdb_sensitivity=mean_pdb_sensitivity,
    mean_pdb_specificity=mean_pdb_specificity,
    mean_pdb_accuracy=mean_pdb_accuracy)
    )
}

write.csv(predictor_summary_scores, file = "predictor_summary_scores.csv")

# correlations from results



#vector
correlations = c()
for (p1 in predictors) {
    s1 = all_results[all_results$predictor==p1,]
    for (p2 in predictors) {
        s2 = all_results[all_results$predictor==p2,]
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

disprot = read.csv("disprot.csv", stringsAsFactors=FALSE)
head(disprot)
#
for (p in predictors) {
    scores = all_results[all_results$predictor==p,c('protein_id', 'position','score')]
    disprot = merge(disprot,scores,by=c('protein_id','position'))
    new_column_names = c(head(colnames(disprot),-1), p)
    colnames(disprot) = new_column_names
}

write.csv(disprot,file='disprot_scores.csv')

# only keep proteins with a PDB, these are more simliar to kinesins
disprot=disprot[disprot$pdb.id!='',]
protein_ids = unique(disprot$protein_id)
# randomly select 95% for training
buckets = sample(2, length(protein_ids), replace=TRUE, prob=c(0.90, 0.10))
protein_buckets = cbind(protein_ids, buckets)
colnames(protein_buckets) = c('protein_id','bucket')
nn2 = merge(disprot,protein_buckets,by=c('protein_id'))
library(class)
# randomly select 95% for training

#bucket = sample(2, nrow(disprot), replace=TRUE, prob=c(0.95, 0.05))

#ones = round(nrow(disprot) * 0.95)
#twos = nrow(disprot) - ones
#bucket = c(rep(1, ones), rep(2, twos))

#training_disordered = nn2[nn2$bucket==1 & nn2$disordered=='Y',]
#training_ordered = nn2[nn2$bucket==1 & nn2$disordered=='N',]
# assume that we have more ordered than disordered, take a subset of ordered to match the size of disordered
#training_ordered2 = training_ordered[sample(nrow(training_ordered), nrow(training_disordered)),]
#training2=rbind(training_disordered, training_ordered2)
#training = training2[,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
#training.labels = training2[,c('disordered')]
training = nn2[nn2$bucket==1,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
training.labels = nn2[nn2$bucket==1,c('disordered')]
testing = nn2[nn2$bucket==2,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
testing.labels = nn2[nn2$bucket==2,c('disordered')]
knn_result = knn(train=training, test=testing, cl=training.labels, k=3)
## write out test results data
testing_df = nn2[nn2$bucket==2,]
testing_df = testing_df[order(testing_df["protein_id"], testing_df["position"]),]
# sort it
knn_results = cbind(testing_df, knn_result)
knn_results_sorted = knn_results[order(knn_results["protein_id"], knn_results["position"]),]
write.csv(knn_results_sorted, file = "knn_results.csv")

true_positive = sum((knn_result=='Y') & (testing.labels =='Y'))
# false negative
false_negative = sum((knn_result=='N') & (testing.labels =='Y'))
# false positive
false_positive = sum((knn_result =='Y') & (testing.labels =='N'))
# true negative
true_negative = sum((knn_result =='N') & (testing.labels =='N'))


knn_sensitivity = true_positive / (true_positive + false_negative)
knn_specificity = true_negative / (true_negative + false_positive)
knn_accuracy = (true_positive + true_negative) / (true_positive + false_positive + false_negative + true_negative)
knn_sensitivity
knn_specificity
knn_accuracy



write.csv(nn2, file = "nn2.csv")

plasmodb_refs = Sys.glob('plasmodb-ref/*csv')
plasmodb_ref_dfs = lapply(plasmodb_refs, function(plasmo) {
    read.csv(plasmo, stringsAsFactors=FALSE)
})
# https://stat.ethz.ch/pipermail/r-help/2010-September/252046.html
plasmodb_ref = do.call("rbind", plasmodb_ref_dfs)
write.csv(plasmodb_ref, file = "plasmodb_ref.csv")
# get predictions
plasmo_prediction_files = Sys.glob('plasmodb/*.csv')
plasmo_prediction_dfs = lapply(plasmo_prediction_files, function(plasmo_prediction_file) {
    read.csv(plasmo_prediction_file, stringsAsFactors=FALSE)
})
# https://stat.ethz.ch/pipermail/r-help/2010-September/252046.html
plasmo_prediction_results = do.call("rbind", plasmo_prediction_dfs)
write.csv(plasmo_prediction_results, file = "plasmo_prediction_results.csv")
plasmodb = plasmodb_ref
for (p in predictors) {
    scores = plasmo_prediction_results[plasmo_prediction_results$predictor==p,c('protein_id', 'position','score')]
    plasmodb = merge(plasmodb,scores,by=c('protein_id','position'))
    new_column_names = c(head(colnames(plasmodb),-1), p)
    colnames(plasmodb) = new_column_names
}
write.csv(plasmodb, file = "plasmodb.csv")

plasmo_knn.training = disprot[,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
plasmo_knn.training_labels = disprot[,c('disordered')]
plasmo_knn.testing = plasmodb[,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]

plasmo_knn.results = knn(train=plasmo_knn.training, test=plasmo_knn.testing, cl=plasmo_knn.training_labels, k=3)

j = cbind(plasmodb, plasmo_knn.results)
j = j[order(j["protein_id"], j["position"]),]
write.csv(j, file = "j.csv")
# with equal
plasmo_knn.training = training
plasmo_knn.training_labels = training.labels
plasmo_knn.testing = plasmodb[,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
plasmo_knn.testing = plasmodb[,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]

plasmo_knn.results = knn(train=plasmo_knn.training, test=plasmo_knn.testing, cl=plasmo_knn.training_labels, k=3)

j = cbind(plasmodb, plasmo_knn.results)
j = j[order(j["protein_id"], j["position"]),]
write.csv(j, file = "j2.csv")

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
# linear regression fit
lm_y = sapply(disprot$disordered, y_score)
lm_y = as.vector(lm_y)
str(disprot$disordered)
str(lm_y)
head(disprot$disordered)
head(lm_y)
disprot_lm = cbind(nn2, lm_y)
# create lm_data
lm_data = disprot_lm[disprot_lm$bucket==1,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short','lm_y')]
# get a fit for lm_y based on prediction scores
lm_fit = lm(lm_y ~ ., data = lm_data)
# prediction
lm_test_data = disprot_lm[disprot_lm$bucket==2,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
lm_test_labels = disprot_lm[disprot_lm$bucket==2,c('disordered')]
lm_predictions = predict(lm_fit, lm_test_data)
lm_predictions_YN = sapply(lm_predictions,to_disorder_YN)
lm_predictions_YN = as.vector(lm_predictions_YN)
lm_results = disprot_lm[disprot_lm$bucket==2,]
lm_results = cbind(lm_results, lm_predictions, lm_predictions_YN)
lm_results = lm_results[order(lm_results["protein_id"], lm_results["position"]),]
write.csv(lm_results, file = "lm_results.csv")

all_results_sorted = cbind(knn_results_sorted, lm_results$lm_predictions_YN)
# write.csv(all_results_sorted, file = "all_results_sorted.csv")
# append lm results

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
lm_sensitivity
lm_specificity
lm_accuracy

library("e1071")
yes_weight=4
num_samples=5000

svm_train.all = nn2[nn2$bucket==1,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short','disordered')]
training_subset = sample(nrow(svm_train.all), num_samples)
svm_train.subset = svm_train.all[training_subset,]
svm_train.data = svm_train.subset[,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
svm_train.labels = as.factor(svm_train.subset[,c('disordered')])

svm_model = svm(x=svm_train.data, y=svm_train.labels, probability = TRUE, class.weights=c('Y'=yes_weight, 'N'=1))
# nrow(svm_train.data.all)

svm_test_all = nn2[nn2$bucket==2,c('protein_id','position','disordered','DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
svm_test_all = svm_test_all[order(svm_test_all["protein_id"], svm_test_all["position"]),]
svm_test_data = svm_test_all[,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
svm_test_labels = svm_test_all[,c('disordered')]
# svm_prediction = predict(svm_model, svm_test_data, decision.values = TRUE, probability = TRUE)
svm_prediction = predict(svm_model, svm_test_data)
svm_prediction_YN = as.vector(svm_prediction)
svm_results = cbind(svm_test_all, svm_prediction_YN)

svm_prediction_filter_lone_Y = apply(svm_results, 1, function(svm_result_row) {
protein_id = svm_result_row['protein_id']
position = as.integer(svm_result_row['position'])
prediction = svm_result_row['svm_prediction_YN']
prediction_previous = svm_results[svm_results$protein_id == protein_id & svm_results$position == (position - 1), c('svm_prediction_YN')]
if (length(prediction_previous) == 0) {
prediction_previous = "_"
}
prediction_next = svm_results[svm_results$protein_id == protein_id & svm_results$position == (position + 1), c('svm_prediction_YN')]
if (length(prediction_next) == 0) {
prediction_next = "_"
}
if (prediction == 'Y' & prediction_previous == 'N' & prediction_next == 'N') {
return ('N')
}
else {
    return (prediction)
}
# check previous prediction, and next prediction
})
svm_results = cbind(svm_results, svm_prediction_filter_lone_Y)

# then filter lone N
svm_prediction_filtered = apply(svm_results, 1, function(svm_result_row) {
    protein_id = svm_result_row['protein_id']
position = as.integer(svm_result_row['position'])
prediction = svm_result_row['svm_prediction_filter_lone_Y']
prediction_previous = svm_results[svm_results$protein_id == protein_id & svm_results$position == (position - 1), c('svm_prediction_filter_lone_Y')]
if (length(prediction_previous) == 0) {
    prediction_previous = "_"
}
prediction_next = svm_results[svm_results$protein_id == protein_id & svm_results$position == (position + 1), c('svm_prediction_filter_lone_Y')]
if (length(prediction_next) == 0) {
    prediction_next = "_"
}
if (prediction == 'N' & prediction_previous == 'Y' & prediction_next == 'Y') {
    return ('Y')
}
else {
    return (prediction)
}
# check previous prediction, and next prediction
})
svm_results = cbind(svm_results, svm_prediction_filtered)

write.csv(svm_results, file = "svm_results.csv")

true_positive = sum((svm_prediction_YN=='Y') & (svm_test_labels =='Y'))
# false negative
false_negative = sum((svm_prediction_YN=='N') & (svm_test_labels =='Y'))
# false positive
false_positive = sum((svm_prediction_YN =='Y') & (svm_test_labels =='N'))
# true negative
true_negative = sum((svm_prediction_YN =='N') & (svm_test_labels =='N'))

svm_sensitivity = true_positive / (true_positive + false_negative)
svm_specificity = true_negative / (true_negative + false_positive)
svm_accuracy = (true_positive + true_negative) / (true_positive + false_positive + false_negative + true_negative)

c('svm_sensitivity'=svm_sensitivity, 'svm_specificity'=svm_specificity, 'svm_accuracy'=svm_accuracy)

## filtered scoring
true_positive = sum((svm_prediction_filtered=='Y') & (svm_test_labels =='Y'))
# false negative
false_negative = sum((svm_prediction_filtered=='N') & (svm_test_labels =='Y'))
# false positive
false_positive = sum((svm_prediction_filtered =='Y') & (svm_test_labels =='N'))
# true negative
true_negative = sum((svm_prediction_filtered =='N') & (svm_test_labels =='N'))

svm_filtered_sensitivity = true_positive / (true_positive + false_negative)
svm_filtered_specificity = true_negative / (true_negative + false_positive)
svm_filtered_accuracy = (true_positive + true_negative) / (true_positive + false_positive + false_negative + true_negative)

c('svm_filtered_sensitivity'=svm_filtered_sensitivity, 'svm_filtered_specificity'=svm_filtered_specificity, 'svm_filtered_accuracy'=svm_filtered_accuracy)



plasmo_svm.test_data = plasmodb[,c('DISEMBL_COILS','DISEMBL_REM465','DISEMBL_HOTLOOPS','DISOPRED','iupred_long','iupred_short')]
plasmo_svm.prediction = predict(svm_model, plasmo_svm.test_data)
plasmo_svm.prediction_YN = as.vector(plasmo_svm.prediction)
plasmo_results = cbind(plasmodb, plasmo_svm.prediction_YN)
plasmo_results = plasmo_results[order(plasmo_results["protein_id"], plasmo_results["position"]),]
write.csv(plasmo_results, file = "plasmo_results.csv")
