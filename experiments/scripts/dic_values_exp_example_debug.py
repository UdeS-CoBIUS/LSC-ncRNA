# This hard coded values values for exeperemenst results for debugging purpose only.
# dictionary structure:{
# 'fixed_length': {
#       dataset_size[50,100,...,350]:{
#           common motif length [2,3,...,20]:{
#               dict of exp results which is
#               {
#                "dataset_size": value must be same as key [50,100,...,350],
#                "min_length": value must be same as key common motif length [2,3,...,20],
#                "max_length": value same as min_length,
#                "num_motifs": value proportional to dataset_size and common motif length => bigger dataset_size and common motif length mean bigger num_motifs,
#                "cpp_execution_time_sec": value proportional to num_motifs,
#                "file_size_gb": value proportional to num_motifs,
#                'EXT': {
#                    "accuracy": proportional to dataset_size, decrease slitly with the increase in dataset_size,
#                    "precision": value,
#                    "recall": value,
#                    "f1_score": value,
#                    "training_time_sec": proportional to dataset_size,
#                    "testing_time_sec": proportional to dataset_size,
#                },
#                'MLP': {
#                    "accuracy": generaly get better with more dataset_size,
#                    "precision": value,
#                    "recall": value,
#                    "f1_score": value,
#                    "training_time_sec": proportional to dataset_size,
#                    "testing_time_sec": proportional to dataset_size,
#                }
#            }
#               }
#           }
#       }
# 'combined_length':{
#       dataset_size[50,100,...,350]:{
#           min_len (fixed at 2), max_len (change from 2 to 20) [(2,2),(2,3),...,(2,20)]:{
#               dict of exp results which is
#               {
#                "dataset_size": value must be same as key [50,100,...,350],
#                "min_length": value must be same as key min_len,
#                "max_length": value must be same as key max_len,
#                "num_motifs": value proportional to dataset_size and common motif length => bigger dataset_size and common motif length mean bigger num_motifs,
#                "cpp_execution_time_sec": value proportional to num_motifs,
#                "file_size_gb": value proportional to num_motifs,
#                'EXT': {
#                    "accuracy": proportional to dataset_size, decrease slitly with the increase in dataset_size,
#                    "precision": value,
#                    "recall": value,
#                    "f1_score": value,
#                    "training_time_sec": proportional to dataset_size,
#                    "testing_time_sec": proportional to dataset_size,
#                },
#                'MLP': {
#                    "accuracy": generaly get better with more dataset_size,
#                    "precision": value,
#                    "recall": value,
#                    "f1_score": value,
#                    "training_time_sec": proportional to dataset_size,
#                    "testing_time_sec": proportional to dataset_size,
#                }
#            }
#               }
#           }
#       }
#}

# todo: accuracy is not good. in EXT it should decrease , but no.

def calculate_metrics(base, dataset_size, motif_length, is_ext=True):
    #dataset_factor = (dataset_size / 350) ** 1.5  # Non-linear scaling for dataset size
    dataset_factor = (dataset_size / 350) **  2.0 # Make scaling more prominent for EXT
    motif_factor = motif_length / 40  # Normalize motif length
    
    if is_ext:
        # EXT accuracy decreases with dataset size
        # accuracy = base - 0.20 * dataset_factor + 0.1 * motif_factor
        accuracy = base - 0.4 * dataset_factor + 0.05 * motif_factor
    else:  # MLP
        # MLP accuracy increases with dataset size
        #accuracy = base + 0.15 * dataset_factor + 0.15 * motif_factor
        accuracy = base + 0.2 * dataset_factor + 0.15 * motif_factor
    
    accuracy = max(0.5, min(0.99, accuracy))  # Clamp between 0.5 and 0.99
    precision = max(0.5, min(0.99, accuracy - 0.02))
    recall = max(0.5, min(0.99, accuracy - 0.01))
    f1_score = max(0.5, min(0.99, (precision + recall) / 2))
    
    return accuracy, precision, recall, f1_score

dict_len_motifs_exp_example_debug = {
    'fixed_length': {
        dataset_size: {
            motif_length: {
                "dataset_size": dataset_size,
                "min_length": motif_length,
                "max_length": motif_length,
                "num_motifs": int(dataset_size * motif_length / 10),
                "cpp_execution_time_sec": dataset_size * motif_length / 100,
                "file_size_gb": dataset_size * motif_length / 10000,
                'EXT': {
                    "accuracy": ext_accuracy,
                    "precision": ext_precision,
                    "recall": ext_recall,
                    "f1_score": ext_f1_score,
                    "training_time_sec": dataset_size * motif_length / 10,
                    "testing_time_sec": dataset_size * motif_length / 40,
                },
                'MLP': {
                    "accuracy": mlp_accuracy,
                    "precision": mlp_precision,
                    "recall": mlp_recall,
                    "f1_score": mlp_f1_score,
                    "training_time_sec": dataset_size * motif_length / 8,
                    "testing_time_sec": dataset_size * motif_length / 35,
                }
            }
            for motif_length in range(2, 21)
            for ext_accuracy, ext_precision, ext_recall, ext_f1_score in [calculate_metrics(0.95, dataset_size, motif_length, True)]
            for mlp_accuracy, mlp_precision, mlp_recall, mlp_f1_score in [calculate_metrics(0.75, dataset_size, motif_length, False)]
        }
        for dataset_size in [50, 100, 150, 200, 250, 300, 350]
    },
    'combined_length': {
        dataset_size: {
            (2, max_len): {
                "dataset_size": dataset_size,
                "min_length": 2,
                "max_length": max_len,
                "num_motifs": int(dataset_size * max_len / 15),
                "cpp_execution_time_sec": dataset_size * max_len / 150,
                "file_size_gb": dataset_size * max_len / 15000,
                'EXT': {
                    "accuracy": ext_accuracy,
                    "precision": ext_precision,
                    "recall": ext_recall,
                    "f1_score": ext_f1_score,
                    "training_time_sec": dataset_size * max_len / 12,
                    "testing_time_sec": dataset_size * max_len / 45,
                },
                'MLP': {
                    "accuracy": mlp_accuracy,
                    "precision": mlp_precision,
                    "recall": mlp_recall,
                    "f1_score": mlp_f1_score,
                    "training_time_sec": dataset_size * max_len / 10,
                    "testing_time_sec": dataset_size * max_len / 40,
                }
            }
            for max_len in range(2, 21)
            for ext_accuracy, ext_precision, ext_recall, ext_f1_score in [calculate_metrics(0.93, dataset_size, max_len, True)]
            for mlp_accuracy, mlp_precision, mlp_recall, mlp_f1_score in [calculate_metrics(0.73, dataset_size, max_len, False)]
        }
        for dataset_size in [50, 100, 150, 200, 250, 300, 350]
    }
}

# Example of how to access the data:
# print(dict_len_motifs_exp_example_debug['fixed_length'][100][5]['EXT']['accuracy'])
# This will print the accuracy for EXT with dataset size 100 and motif length 5 in the fixed_length scenario
"""

1. For EXT:
   - Old: `"accuracy": min(0.95, 0.85 + dataset_size / 1000)`
   - New: `"accuracy": max(0.85, 0.95 - dataset_size / 3500)`
   
   This change makes the accuracy start high (0.95) and decrease slightly as the dataset size increases, with a minimum of 0.85. This reflects the description that EXT's accuracy is proportional to dataset size but decreases slightly with larger dataset sizes.

2. For MLP:
   - Kept as: `"accuracy": min(0.93, 0.80 + dataset_size / 1000)`
   
   This formula remains unchanged as it already correctly implements the behavior where MLP's accuracy improves with larger dataset sizes, starting at 0.80 and increasing up to a maximum of 0.93.

With these changes:
- EXT's accuracy will start higher but decrease slightly as the dataset size increases.
- MLP's accuracy will start lower but increase as the dataset size increases.
- For very large dataset sizes, MLP's accuracy may surpass EXT's, which aligns with the idea that MLP generally improves with more data.

These modifications ensure that the dictionary now correctly represents the described behavior for both EXT and MLP accuracies. The same logic has been applied to precision, recall, and f1_score to maintain consistency.

Is this modification more in line with what you were expecting? Would you like any further adjustments or explanations?
"""