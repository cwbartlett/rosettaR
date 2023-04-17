library(lavaan)
library(ggplot2)
library(dplyr)

row_set=seq(500,10000,
            length.out=5)

n_row_set = length(row_set)

df = NULL

for(s in row_set){

  n_reps = 50
  results = NULL
  for(r in 1:n_reps){
    d <- sim(n_rows = s)

    # run rosetta
    d_rosetta <- rosetta(
      d = d$missing,
      factor_structure = list(
        a = c("a_1", "a_2", "a_3"),
        b = c("b_1", "b_2", "b_3"),
        c = c("c_1", "c_2", "c_3")
      )
    )

    obs_cov = get_obs_cov(rosetta_bind(d$missing))[c("a_1",
                                                     "a_2",
                                                     "a_3",
                                                     "b_1",
                                                     "b_2",
                                                     "b_3",
                                                     "c_1",
                                                     "c_2",
                                                     "c_3"),
                                                   c("a_1",
                                                      "a_2",
                                                      "a_3",
                                                      "b_1",
                                                      "b_2",
                                                      "b_3",
                                                      "c_1",
                                                      "c_2",
                                                      "c_3")]

    mod_obj  = d_rosetta |> attr( "unconstrained_fit_lavaan_object")
    pred_cov <- matrix(fitted.values(mod_obj) |> unlist(),nrow=nrow(obs_cov))
    error=(pred_cov-obs_cov)^2
    error = sum(error[upper.tri(error,diag = TRUE)])
    results= rbind(results,c(error))
  }
  df=rbind(df,
           c(s,mean((results)),var((results)),"matrix_error"))
}

colnames(df) <- c("n_rows","mse","vse","par_name")
df = df |> as.data.frame()


df$n_rows <- as.numeric(df$n_rows)
df$mse <- as.numeric(df$mse)
df$vse <- as.numeric(df$vse)

require(ggplot2)
ggplot(df) + geom_line(aes(x=n_rows,y = mse ,color=par_name)) +
  geom_ribbon(aes(x=n_rows,ymax = mse + sqrt(vse/n_reps),
                  ymin=  mse - sqrt(vse/n_reps),fill=par_name),alpha=.1)

