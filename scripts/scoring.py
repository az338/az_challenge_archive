import subprocess

# path to modified R scoring script for challenge 2
scoring_script=
# path to predicted combinations file with prediction scores (matrix)
pred_cc = 
# path to combination priorities file with confidence score (matrix)
pred_confidence =
# path to observed combinations file with synergy scores  (melted data frame with one score per row)
leaderboard_cc = 

# runs R challenge 2 scoring script with above files
subprocess.call(['Rscript',scoring_script+'',pred_cc,pred_confidence,leaderboard_cc])
