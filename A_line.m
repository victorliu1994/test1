function [SNR, noise_power, computation_time, RSS, recovery, original] = A_line(B_scan_path, A_line_index, sampling_type, sampling_rate, bases, noise_mask, object_function_parameter)
    tic;
    img = double(imread(B_scan_path));
    %% Noise&SNR compute part
    noise_est_area = double(img(noise_mask(1):noise_mask(2),noise_mask(3):noise_mask(4)));
    noise_power = 10*log10(var(noise_est_area(:)));
    med_filter_res =  medfilt2(img,[5,5]);
    SNR = 10*log10(var(med_filter_res(:))/var(noise_est_area(:)));
    %% Compress Sensing part
    original_signal = img(:,A_line_index);
    basis_transform_matrix = zeros(length(original_signal),length(original_signal));
    if(strcmp(bases,'default') || strcmp(bases,'wave'))
        basis_transform_matrix = DWT_test(length(original_signal));
    elseif(strcmp(bases,'fourier'))
        basis_transform_matrix = dctmtx(length(original_signal));
    else
         error("Wrong bases: No such bases for current version, may extend further");
    end
    
    len = length(object_function_parameter);
    if(len < 1)
        object_function_coef(1) = 1;
    elseif(len == 1)
        object_function_coef(1) = object_function_parameter(1);
    elseif(len == 2)
        object_function_coef(1) = object_function_parameter(1);
        object_function_coef(2) = object_function_parameter(2);
    elseif(len == 3)
        object_function_coef(1) = object_function_parameter(1);
        object_function_coef(2) = object_function_parameter(2);
        object_function_coef(3) = object_function_parameter(3);
    elseif(len == 4)
        object_function_coef = object_function_parameter;
    else
        error("Wrong object_function_parameter: too much coefficients for object function");
    end
    n = length(original_signal);
    sampling_matrix = zeros(sampling_rate, n);
    if(strcmp(sampling_type,'default') || strcmp(sampling_type,'gaussian random'))
        sampling_matrix = randn(sampling_rate, n);
    elseif(strcmp(sampling_type, 'uniform'))
        for i = 1:n
            for j = 1:sampling_rate
                if (i/j == n/sampline_rate)
                    sampling_matrix(j, i) = 1; 
                end
            end
        end
    else
        error("Wrong sampling_type: no such sampling type for current version, may extend further");
    end
    
    compress_signal = sampling_matrix*basis_transform_matrix*original_signal;
    %% Recovery and calculate residual sum of squares
    cvx_begin
        variable x(n)
        minimize (object_function_coef(1)*norm(x,1)+object_function_coef(2)*norm(x,2)+object_function_coef(3)*norm_nuc(x))
        subject to
            sampling_matrix*x == compress_signal
    cvx_end
    original = original_signal;
    recovery = inv(basis_transform_matrix)*x;
    RSS = sum((original-recovery).*(original-recovery));
    toc;
    computation_time = toc;
end
