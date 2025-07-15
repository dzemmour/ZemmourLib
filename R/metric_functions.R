#' Calculate True Positives, True Negatives, False Positives, and False Negatives
#'
#' This helper function calculates the counts of TP, TN, FP, and FN given two
#' binary vectors (e.g., one from a gating strategy and one from a true cluster).
#'
#' @param x A logical vector (e.g., `colnames(so) %in% gate_list[["cl2"]]`).
#' @param y A logical vector (e.g., `so$is_cl2`).
#' @return A named numeric vector with counts for tp, tn, fp, fn.
#' @examples
#' tp_tn_fp_fn_example <- function() {
#'   set.seed(123)
#'   true_status <- sample(c(TRUE, FALSE), 100, replace = TRUE, prob = c(0.3, 0.7))
#'   gate_result <- sample(c(TRUE, FALSE), 100, replace = TRUE, prob = c(0.25, 0.75))
#'   TpTnFpFn(gate_result, true_status)
#' }
#' tp_tn_fp_fn_example()
#' @export
TpTnFpFn = function(x, y) {
    tmp = table(x,y)
    # Ensure all combinations are present, filling with 0 if not
    tmp_matrix <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("FALSE", "TRUE"), c("FALSE", "TRUE")))
    tmp_matrix[rownames(tmp), colnames(tmp)] <- tmp

    tp = tmp_matrix["TRUE","TRUE"]
    tn = tmp_matrix["FALSE","FALSE"]
    fp = tmp_matrix["TRUE","FALSE"]
    fn = tmp_matrix["FALSE","TRUE"]
    return(c(tp = tp, tn = tn, fp=fp, fn=fn))
}

#' Calculate Sensitivity, Specificity, PPV, and NPV
#'
#' This function calculates sensitivity (Se), specificity (Sp), positive predictive value (PPV),
#' and negative predictive value (NPV) from true positive, true negative, false positive, and false negative counts.
#'
#' @param tp Numeric, true positives count.
#' @param tn Numeric, true negatives count.
#' @param fp Numeric, false positives count.
#' @param fn Numeric, false negatives count.
#' @return A named numeric vector with se, sp, ppv, npv.
#' @examples
#' counts <- TpTnFpFn(x = c(T,T,F,F,T), y = c(T,F,T,F,T))
#' SeSpPPVNPV(tp = counts["tp"], tn = counts["tn"], fp = counts["fp"], fn = counts["fn"])
#' @export
SeSpPPVNPV = function(tp = NULL, tn = NULL, fp = NULL, fn = NULL) {
    if (is.null(tp) || is.null(tn) || is.null(fp) || is.null(fn)) {
        stop("All arguments (tp, tn, fp, fn) must be provided.")
    }
    # Handle division by zero gracefully, return NA if denominator is zero
    se = ifelse((tp + fn) > 0, tp / (tp+fn), NA)
    sp = ifelse((fp + tn) > 0, tn / (fp+tn), NA)
    ppv = ifelse((tp + fp) > 0, tp / (tp+fp), NA)
    npv = ifelse((tn + fn) > 0, tn / (tn+fn), NA)
    res = c(se, sp, ppv, npv)
    names(res) =  c("se", "sp", "ppv", "npv")
    return(res)
}
