% This script highlights the importance of accounting for multiple
% hypothesis testing when performing statistical analyses. It generates
% purely random survival and gene expression data for a number of patients,
% and performs survival analysis by calculating the Cox proportional 
% hazards p-value for each gene. Histograms of a) the unadjusted p-values,
% and b) the Benjamini-Hochberg p-values are displayed, which demonstrates
% that a failure to account for multiple hypothesis testing can enable type
% I errors.

rng(42);

%% Generate Survival Data

% Generate random survival data for 500 cancer patients over a span of 2000
% days. Assume that on each day, there is a 1/2000 risk of death, and a 
% 1/2000 chance of follow-up being lost.

NUM_PTS = 500;
NUM_DAYS = 2000;
DEATH_RISK = 1/2000;
CENSOR_RISK = 1/2000;

days = zeros(NUM_PTS, 1);
outcomes = zeros(NUM_PTS, 1);

for pt = 1:NUM_PTS
    [days(pt), outcomes(pt)] = generate_pt_survival_data(NUM_DAYS, DEATH_RISK, CENSOR_RISK); 
end

surv_data = sortrows(table(days, outcomes, 'VariableNames', {'Days', 'IsDeceased'}));

clear NUM_DAYS DEATH_RISK CENSOR_RISK days outcomes pt

%% Generate Expression Data

% Generate random expression data for each patient for 100 different
% 'genes'.
NUM_GENES = 100;
MEAN_EXPS = repmat(101:(100 + NUM_GENES), NUM_PTS, 1);
MEAN_SDS = 0.2 * MEAN_EXPS;

% Construct an expression table with more suitable variable names.
expr_data = array2table(MEAN_EXPS + MEAN_SDS .* randn(NUM_PTS, NUM_GENES));
expr_data.Properties.VariableNames = regexprep(expr_data.Properties.VariableNames, 'Var', 'Gene');

clear MEAN_EXPS MEAN_SDS

%% Perform univariate Cox proportional hazards regression for survival analysis

cox_pvals = zeros(NUM_GENES, 1);

for gene = 1:NUM_GENES
    [~, ~, ~, stats] = coxphfit(expr_data{:, gene}, surv_data.Days, 'Censoring', ~surv_data.IsDeceased);
    cox_pvals(gene) = stats.p;
end

figure
histogram(cox_pvals, 0:0.05:1, 'Normalization', 'probability');
ylabel('Proportion of Genes');
xlabel('Univariate Cox PH p-value');

fprintf('Without adjustment, %i of the %i genes have p-values that are significant - even though the data were randomly generated!\n', sum(cox_pvals < 0.05), NUM_GENES);

clear gene stats

%% Adjust the p-values for multiple hypothesis testing

adj_pvals = sort(cox_pvals);
adj_pvals = adj_pvals * NUM_GENES ./ (1:NUM_GENES)';

for val = NUM_GENES:-1:2
    if adj_pvals(val) < adj_pvals(val - 1)
        adj_pvals(val - 1) = adj_pvals(val);
    end
end

adj_pvals(adj_pvals > 1) = 1;

figure
histogram(adj_pvals, 0:0.05:1, 'Normalization', 'probability');
ylabel('Proportion of Genes');
xlabel('Benjamini-Hochberg Adjusted Univariate Cox PH p-value');

fprintf('After Benjamini-Hochberg adjustment for multiple hypothesis testing:\n%i of the %i genes have p-values that are significant.\n', sum(adj_pvals < 0.05), NUM_GENES);

clear val

%% Functions

function [day, outcome] = generate_pt_survival_data(total_days, risk_of_death, risk_of_censoring)

    raw_surv_data = rand(total_days, 1);

    day_of_first_event = find(raw_surv_data <= (risk_of_death + risk_of_censoring), 1);

    if(~isempty(day_of_first_event))
        day = day_of_first_event;
        outcome = raw_surv_data(day_of_first_event) <= risk_of_death;
    else
        % If the patient didn't have an event, assume that follow-up is not
        % maintained beyond the final day.
        day = total_days;
        outcome = 0;
    end
end