function [features, valid] = features_detection(Im)   

    if ndims(Im) == 3
        I = rgb2gray(Im);
    else
        I = Im;
    end
    
%     %% method 1: Harris
%     corners = detectHarrisFeatures(I);
%     [features, valid] = extractFeatures(I, corners);
    
    %% method 2: SURF
    points = detectSURFFeatures(I);
    [features, valid] = extractFeatures(I, points);

%     %% method 3: MSER
%     regions = detectMSERFeatures(I);
%     [features, valid] = extractFeatures(I,regions,'Upright',true);
    
end