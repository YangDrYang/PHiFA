function compareMeasurements(vec_x, vec_1, vec_2, y_label, x_label)

if length(vec_x)~=length(vec_2) || length(vec_x)~=length(vec_1)
    fprintf('Can not compare measurements. Length are unequal.\n');
    return;
end

n_vec = min([size(vec_1,1), size(vec_2,1)]);

figure;
subplot(2,1,1);
for i = 1:n_vec
    h = plot(vec_x, vec_1(i,:));
    set(h, 'DisplayName', ['a' num2str(i)]);
    hold on;
    h = plot(vec_x, vec_2(i,:), '--');
    set(h, 'DisplayName', ['b' num2str(i)]);
end
legend show;
xlabel(x_label);
ylabel(y_label);
grid;
hold off;

subplot(2,1,2);
for i = 1:n_vec
    h = plot(vec_x, vec_2(i,:)-vec_1(i,:));
    set(h, 'DisplayName', num2str(i));
    hold on;
end
legend show;
xlabel(x_label);
ylabel(['Difference ' y_label]);
grid;
hold off;

end