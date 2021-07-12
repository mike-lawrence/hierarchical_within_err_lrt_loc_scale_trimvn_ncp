#preamble (options, installs, imports & custom functions) ----

options(warn=1) #really should be default in R
`%!in%` = Negate(`%in%`) #should be in base R!

# specify the packages used:
required_packages = c(
	'github.com/mike-lawrence/aria/aria' # for Stan stuff
	, 'tidyverse' #for all that is good and holy
)

#load the helper functions:
# load the helper functions:
for(file in fs::dir_ls('r')){
	cat('Loading function: ',fs::path_ext_remove(fs::path_file(file)),'()\n',sep='')
	source(file)
}

#install any required packages not already present
install_if_missing(required_packages)

#load packages
library(tidyverse)
library(aria)

#import & toss too-fast
#  Demo note: just data I had on hand for this
(
	readRDS('rds/00_imported_and_cleaned.rds')
	%>% filter(!str_detect(experiment,'ANTI_Poder_Legitimidad')) #no error data in this experiment
	%>% filter(flanker!='neutral') # for demo
	%>% filter(rt>150)
	%>% mutate(lrt=log(rt))
	%>% rename(id=participant)
	%>% filter(id %in% unique(id)[1:30]) #just a few Ss for this demo
) ->
	a

# quick viz:
(
	a
	# %>% filter(lrt<6)
	%>% ggplot()
	+ facet_wrap(
		~error
		, ncol = 1
		, scales = 'free_y'
	)
	+ geom_histogram(
		mapping = aes(
			x = lrt
		)
	)
)


#set factor level order for interpretability of helmert contrasts
(
	a
	%>% mutate(
		alertness = factor(alertness,levels=c('tone','no tone'))
		, orienting = factor(orienting,levels=c('valid','invalid','none'))
		, flanker = factor(flanker,levels=c('congruent','incongruent'))
	)
) -> a


# Compute inputs to model ----

(
	a
	#add the contrast matrix columns
	%>% mutate(
		contrasts = get_contrast_matrix_rows_as_list(
			data = .
			, formula = ~ flanker # Demo note: we could do `~alertness*orienting*flanker`, but that'll take a long time to sample
			, contrast_kind = halfhelmert_contrasts
		)
	)
) -> a

#first get the unique contrasts (will serve as input to Stan)
(
	a
	%>% distinct(contrasts)
	%>% arrange(contrasts) #important
) ->
	unique_contrasts


#now cross with unique subjects
(
	unique_contrasts
	%>% expand_grid(
		tibble(id = unique(a$id))
	)
	%>% arrange(contrasts,id) #important
	# add a column containing the row index
	%>% mutate(
		row = 1:n()
	)
	#join with dat_summary
	%>% right_join(
		a
		, by = c('contrasts','id')
	)
) -> a

#now get the unique contrast matrix
(
	unique_contrasts
	%>% unnest(contrasts)
	%>% as.matrix()
) -> uw


(
	a
	%>% group_by(row)
	%>% summarise(
		cell_count_total = n()
		, cell_count_errors = sum(error)
	)
	%>% arrange(row)
) -> err_stats



# package for stan & sample ----

data_for_stan = lst( #lst permits later entries to refer to earlier entries

	####
	# Entries we need to specify ourselves
	####

	# uw: unique within predictor matrix
	w = uw

	# lrt: lrt observations
	, lrt = a$lrt[!a$error]

	# lrt_row: lrt row indicator
	, lrt_row = a$row[!a$error]

	# err_num_total: number of observations in each subject/condition
	, err_cell_count_total = err_stats$cell_count_total

	# err_num_successes: number of 1's in each subject/condition
	, err_cell_count_errors = err_stats$cell_count_errors

	# err_row: err row indicator
	, err_row = err_stats$row

	####
	# Entries computable from the above
	####

	, k = ncol(uw)
	, n = length(unique(a$id))
	, m_lrt = length(lrt_row)
	, m_err = length(err_row)
	, is_intercept = guess_if_intercept(uw)

)

#quick view
glimpse(data_for_stan)

stan_code_path = 'stan/err_lrt_loc_scale_trimvn_ncp.stan'

#take a look at the model
file.edit(stan_code_path)
#when you click save, aria will syntax check and compile

#alternatively:
aria:::check_syntax_and_maybe_compile(stan_code_path)

# compose
aria::compose(
	data = data_for_stan
	, code_path = stan_code_path
	, out_path = 'nc/samples.nc'
)

# Check & viz posterior ----

post = aria::coda('nc/samples.nc')

# Check treedepth, divergences, & rebfmi
(
	post$draws(group='sample_stats')
	%>% posterior::as_draws_df()
	%>% group_by(.chain)
	%>% summarise(
		max_treedepth = max(treedepth)
		, num_divergent = sum(divergent)
		, rebfmi = var(energy)/(sum(diff(energy)^2)/n()) #n.b. reciprocal of typical EBFMI, so bigger=bad, like rhat
	)
)

# gather summary for core parameters (inc. r̂ & ess)
(
	post$draws(groups='parameters')
	%>% posterior::summarise_draws(.cores=parallel::detectCores())
) ->
	par_summary

# show the ranges of r̂/ess's
(
	par_summary
	%>% select(rhat,contains('ess'))
	%>% summary()
)

#View those with suspect r̂
(
	par_summary
	%>% filter(rhat>1.01)
	%>% View()
)

# Viz means
(
	post$draws(variables=c('err_m','mrt_m','srt_m'))
	%>% posterior::as_draws_array()
	%>% bayesplot::mcmc_intervals()
)

# Viz condition expectations for mean log-RT
(
	post$draws(variables='lrt_emean')
	%>% posterior::as_draws_array()
	%>% bayesplot::mcmc_intervals()
)
