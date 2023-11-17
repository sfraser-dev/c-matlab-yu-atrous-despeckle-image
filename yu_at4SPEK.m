% cleaned_image=yu_at4SPEK(image,th_seed,P,sstep,gamma,filterType,th0Str);
%
% image=noisy input image.
% th_seed=initial threshold (seed).
% P=percentage of noise to be removed from input image.
% sstep=step size.
% gamma=alteration of noise estimate.
% filterType='ST', 'HT' or 'SHT'; hard, soft and soft-hard combination thresholding.
% t0Str=chop off or filter the finest wavelet co-effs ('th0+' and 'th0-' respectively).
%
% speckle corrupted images only.

function cleaned_image=yu_at4SPEK(seedVal,image,th_seed,P,sstep,gamma,filterType,th0Str);
I=double(image); 
[no_rows,no_cols]=size(I);

% Log transform.
I=log(I+1);

% Computing 4 level wavelet decomposition.
s0=hole_wavelet2d(I,0);
s1=hole_wavelet2d(s0,1);
s2=hole_wavelet2d(s1,2);
s3=hole_wavelet2d(s2,3);
w0=I-s0;
w1=s0-s1;
w2=s1-s2;
w3=s2-s3;

% Get a NOISE ONLY image.
noise_only=get_noise_only_image(seedVal,w0,gamma);
%randn('seed',seedVal); noise_only=randn(size(w0))*std(w0(:))*gamma;

% Passing the NOISE ONLY image through the wavelet transform.
ns0=hole_wavelet2d(noise_only,0);
ns1=hole_wavelet2d(ns0,1);
ns2=hole_wavelet2d(ns1,2);
ns3=hole_wavelet2d(ns2,3);
nw0=noise_only-ns0;
nw1=ns0-ns1;
nw2=ns1-ns2;
nw3=ns2-ns3;

% Perform the iterative filtering upon each wavelet level.
% Using soft thresholding.
if strcmp(filterType,'ST')
	fprintf('\n----------level 0----------\n');
	if strcmp(th0Str,'th0+')==1
		[C0,NC0,th0]=ST_iterate(w0,nw0,th_seed,sstep,P);
	elseif strcmp(th0Str,'th0-')==1
		fprintf('skipping...will chop off later\n');
	end
	fprintf('\n----------level 1----------\n');
	[C1,NC1,th1]=ST_iterate(w1,nw1,th_seed,sstep,P);
	fprintf('th1=%f\n',th1);
	fprintf('\n----------level 2----------\n');
	[C2,NC2,th2]=ST_iterate(w2,nw2,th_seed,sstep,P);
	fprintf('th2=%f\n',th2);
	fprintf('\n----------level 3----------\n');
	[C3,NC3,th3]=ST_iterate(w3,nw3,th_seed,sstep,P);
	fprintf('th3=%f\n',th3);
	% Inverse wavelet transform. 
	if strcmp(th0Str,'th0+')==1
		cleaned_image=s3+C3+C2+C1+C0;
	elseif strcmp(th0Str,'th0-')==1 
		cleaned_image=s3+C3+C2+C1;
	end	
% Using hard thresholding.
elseif strcmp(filterType,'HT')
	fprintf('\n----------level 0----------\n');
	if strcmp(th0Str,'th0+')==1
		[C0,NC0,th0]=HT_iterate(w0,nw0,th_seed,sstep,P);
	elseif strcmp(th0Str,'th0-')==1
		fprintf('skipping...will chop off later\n');
	end
	fprintf('\n----------level 1----------\n');
	[C1,NC1,th1]=HT_iterate(w1,nw1,th_seed,sstep,P);
	fprintf('th1=%f\n',th1);
	fprintf('\n----------level 2----------\n');
	[C2,NC2,th2]=HT_iterate(w2,nw2,th_seed,sstep,P);
	fprintf('th2=%f\n',th2);
	fprintf('\n----------level 3----------\n');
	[C3,NC3,th3]=HT_iterate(w3,nw3,th_seed,sstep,P);
	fprintf('th3=%f\n',th3);
	% Inverse wavelet transform. 
	if strcmp(th0Str,'th0+')==1
		cleaned_image=s3+C3+C2+C1+C0;
	elseif strcmp(th0Str,'th0-')==1 
		cleaned_image=s3+C3+C2+C1;
	end	
% Using soft and hard thresholding.
elseif strcmp(filterType,'SHT')
	fprintf('\n----------level 0----------\n');
	if strcmp(th0Str,'th0+')==1
		[C0,NC0,th0]=ST_iterate(w0,nw0,th_seed,sstep,P);
	elseif strcmp(th0Str,'th0-')==1
		fprintf('skipping...will chop off later\n');
	end
	fprintf('\n----------level 1----------\n');
	[C1,NC1,th1]=HT_iterate(w1,nw1,th_seed,sstep,P);
	fprintf('th1=%f\n',th1);
	fprintf('\n----------level 2----------\n');
	[C2,NC2,th2]=HT_iterate(w2,nw2,th_seed,sstep,P);
	fprintf('th2=%f\n',th2);
	fprintf('\n----------level 3----------\n');
	[C3,NC3,th3]=HT_iterate(w3,nw3,th_seed,sstep,P);
	fprintf('th3=%f\n',th3);
	% Inverse wavelet transform.
	if strcmp(th0Str,'th0+')==1
		cleaned_image=s3+C3+C2+C1+C0;
	elseif strcmp(th0Str,'th0-')==1 
		cleaned_image=s3+C3+C2+C1;
	end	
else
	fprintf('ERROR! please indicate hard or soft thresholding.\n');
end

% Inverse log transform.
cleaned_image=exp(cleaned_image)-1;

%imshow([uint8(image) uint8(cleaned_image)]);title('original : cleaned');

%---------------------------------------%
function noisy=get_noise_only_image(seedVal,w,gamma);
% Estimate the std of the noise.
noise_std=w.*(abs(w) < (std(w(:))*3)); 	% Sigma-3-clipping.
noise_std=std(noise_std(:));
%fprintf('\nnoise (std) in input image: %f\n',noise_std);
randn('seed',seedVal);
noisy=randn(size(w))*noise_std*gamma;	% Constructed noise image.

%---------------------------------------%
function [C,NC,th]=ST_iterate(w,nw,th,sstep,P);
% Main component of the algorithm.
fprintf('initial threshold: %f\n',th);
tolerance=std(nw(:))-((std(nw(:))/100)*P);
fprintf('target std for noisy co-effs: %f\n',std(nw(:))-tolerance);
while 1,
	% Getting the noisy co-effs. 
	tw=softThreshold(w,th);
	NC=w-tw;	% Removed noisy co-effs.
	
	% Break out of while loop if inside tolerance.	
	fprintf('std of noisy co-effs: %f\n',std(NC(:)));
	delta=std(nw(:))-std(NC(:));
	if delta <= tolerance
		C=softThreshold(w,th);	% Cleaned co-effs.
		break;
	end
	% Apply a decaying increment to the threshold.
	th=th+sstep*delta;
	fprintf('incremented threshold=%f\n',th);
end

%---------------------------------------%
function [C,NC,th]=HT_iterate(w,nw,th,sstep,P);
% Main component of the algorithm.
fprintf('initial threshold: %f\n',th);
tolerance=std(nw(:))-((std(nw(:))/100)*P);
fprintf('target std for noisy co-effs: %f\n',std(nw(:))-tolerance);
while 1,
	% Getting the noisy co-effs. 
	NC=w.*(abs(w)<=th); 	% Removed noisy co-effs.
	
	% Break out of while loop if inside tolerance.	
	fprintf('std of noisy co-effs: %f\n',std(NC(:)));
	delta=std(nw(:))-std(NC(:));
	if delta <= tolerance
		C=hardThreshold(w,th);	% Cleaned co-effs.
		break;
	end
	% Apply a decaying increment to the threshold.
	th=th+sstep*delta;
	fprintf('incremented threshold=%f\n',th);
end

%---------------------------------------%
% Copied exactly from WaveLab's SoftThresh function
function x = softThreshold(y,t)
res = (abs(y) - t);
res = (res + abs(res))/2;
x   = sign(y).*res;

%---------------------------------------%
% Copied exactly from WaveLab's HardThresh function
function x = hardThreshold(y,t)
x   = y .* (abs(y) > t);
