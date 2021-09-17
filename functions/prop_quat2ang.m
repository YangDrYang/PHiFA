function [prop_new] = prop_quat2ang(prop_old)

prop_new = zeros(13, size(prop_old,2));
prop_new(1:7,:) = prop_old(1:7,:);
prop_new(8:10,:) = (quat2eul(quaternion(norm_quat(prop_old(8:11,:))'))');
prop_new(11:end,:) = prop_old(12:end,:);

end