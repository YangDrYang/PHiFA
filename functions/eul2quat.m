function q = eul2quat( eul, varargin )
%EUL2QUAT Convert Euler angles to quaternion


% Pre-allocate output
q = zeros(4, size(eul,1), 'like', eul);

% Compute sines and cosines of half angles
c = cos(eul'/2);
s = sin(eul'/2);

% The parsed sequence will be in all upper-case letters and validated

        % Construct quaternion
        q = [c(:,1).*c(:,2).*c(:,3) - s(:,1).*s(:,2).*s(:,3); ...
            s(:,1).*c(:,2).*c(:,3) + c(:,1).*s(:,2).*s(:,3); ...
            -s(:,1).*c(:,2).*s(:,3) + c(:,1).*s(:,2).*c(:,3); ...
            c(:,1).*c(:,2).*s(:,3) + s(:,1).*s(:,2).*c(:,3)];


end

