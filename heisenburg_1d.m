%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               DMRG for the  1D Heisenberg Model
%               After S.R.White, Phys. Rev. Lett. 69, 2863 (1992)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Number of states kept m
m=20;
% Number of iterations.  Final lattice size is 2*Niter + 2
NIter = 1000;
% exact energy, for comparison
ExactEnergy = -log(2) + 0.25;
 
fprintf('Iter\tEnergy\t\tBondEnergy\tEnergyError\tTrunc\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Intialize local operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
I= eye(2);
Sz = [1/2 0;0 -1/2];
Sp = [0 0;1 0];
Sm = [0 1 ;0 0];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Initial blocks
%               We assume reflection symmetry so we only need 1 block
%               The operator acts on the inner-most site of the block
%               +---------+    +---------+
%               |        *|    |*        |
%               +---------+    +---------+
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
BlockSz = Sz;
BlockSp = Sp;
BlockSm = Sm;
BlockI  = I;
BlockH  = zeros(2);
Energy = -0.75;   % energy of 2-site lattice
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Begin main iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for l = 1:NIter
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           Get the 2m-dimensional operators for the block + site
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    BlockH = kron(BlockH, I) + kron(BlockSz, Sz) ...
        + 0.5 * ( kron(BlockSp, Sm) + kron(BlockSm, Sp) );
    BlockSz = kron(BlockI, Sz);
    BlockSp = kron(BlockI, Sp);
    BlockSm = kron(BlockI, Sm);
    BlockI  = kron(BlockI, I);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   HAMILTONIAN MATRIX forsuperblock
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    H_super = kron(BlockH, BlockI) + kron(BlockI, BlockH) ...
        + kron(BlockSz, BlockSz) ...
        + 0.5 * ( kron(BlockSp, BlockSm) + kron(BlockSm, BlockSp) );
 
    H_super = 0.5 * (H_super + H_super');  % ensure H is symmetric
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Diagonalizing the Hamiltonian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    LastEnergy = Energy;
    opts.disp = 0;   % disable diagnostic information in eigs
    opts.issym = 1;
    opts.real = 1;
    [Psi Energy] = eigs(H_super,1,'SA', opts);
    EnergyPerBond = (Energy - LastEnergy) / 2;
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Form the reduced density matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    [nr,nc]=size(Psi);
    Dim = sqrt(nr);
    PsiMatrix = reshape(Psi, Dim, Dim);
    Rho = PsiMatrix * PsiMatrix';
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Diagonalize the density matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    [V D] = eig(Rho);
    [D Index] = sort(diag(D), 'descend');  % sort eigenvalues descending
    V = V(:,Index);                     % sort eigenvectors the same way
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Construct the truncation operator
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    NKeep = min(size(D, 1), m);
    T = V(:, 1:NKeep);
    TruncationError = 1 - sum(D(1:NKeep));
 
    fprintf('%d\t%f\t%f\t%f\t%f\n',l, Energy, EnergyPerBond, ...
        ExactEnergy-EnergyPerBond, TruncationError);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Transform the block operators into the truncated basis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    BlockH  = T'*BlockH*T;
    BlockSz = T'*BlockSz*T;
    BlockSp = T'*BlockSp*T;
    BlockSm = T'*BlockSm*T;
    BlockI = T'*BlockI*T;
 
end