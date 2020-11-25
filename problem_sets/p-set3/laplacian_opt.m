function Lap = laplacian_opt(u,x,y)
    % optimized computation of the Laplacian using Chebyshev grid
    % Instead of doing a for loop we use broadcasting to fully vectorize
    % the operations
    N = size(x,2)-1;

    % we allocate the zeros
    uxx = zeros(N+1,N+1); uyy = zeros(N+1,N+1);

    % we create these indices for simplicity
    ii = 2:N;

    % 2nd derivs wrt x in each row
    V = [u fliplr(u(:,ii))];       % we flip it as a tensor
    U = real(fft(V, [], 2)); % FFT across the rows
    % diff wrt theta
    W1 = real(ifft(bsxfun(@times, 1i*[0:N-1 0 1-N:-1], U), [], 2)); 
    % diff^2 wrt theta
    W2 = real(ifft(bsxfun(@times,-[0:N 1-N:-1].^2, U), [], 2)); 
    % computing the chain rule
    uxx(ii,ii) =  bsxfun(@times,W2(ii,ii),1./(1-x(ii).^2)) ...
                - bsxfun(@times,W1(ii,ii),x(ii)./(1-x(ii).^2).^(3/2));


    % 2nd derivs wrt y in each column
    V = [u; flipud(u(ii,:))];
    U = real(fft(V,[],1)); % FFT across the columns
    % diff wrt theta
    W1 = real(ifft(bsxfun(@times, 1i*[0:N-1 0 1-N:-1].', U), [], 1)); 
    % diff^2 wrt theta
    W2 = real(ifft(bsxfun(@times,(-[0:N 1-N:-1].^2).', U), [], 1));  
    % computing the chain rule
    uyy(ii,ii) =  bsxfun(@times,W2(ii,ii),1./(1-y(ii).^2)) ...
                - bsxfun(@times,W1(ii,ii),y(ii)./(1-y(ii).^2).^(3/2));

    Lap = uxx + uyy;

end