function sim_out_new = sort_sim_out_radau(sim_out)

tepoch = zeros(length(sim_out),1);
for i = 1:length(sim_out)
    tepoch(i) = sim_out(i).t;
end
[~, ind] = sort(tepoch);

sim_out_new(1:length(sim_out),1) = clSimulationOutput();
tepoch_new = zeros(length(sim_out),1);
for i = 1:length(sim_out)
    sim_out_new(i) = sim_out(ind(i));
    tepoch_new(i) = sim_out(ind(i)).t;
end

[~,ind,~] = unique(tepoch_new,'last');
sim_out_new = sim_out_new(ind);







