function run_colpitts_hopf
%RUN_COLPITTS_HOPF  Continuation workflow for the Colpitts oscillator in MatCont.
%
%  This driver locates the Hopf bifurcation on the equilibrium branch while
%  varying g, evaluates the first Lyapunov coefficient, initializes a limit-cycle
%  continuation from the Hopf point, and inspects Floquet multipliers to infer
%  stability. The procedure is repeated over a grid of Q values in [0.5, 3].

init;
odefile = @colpitts;
qValues = 0.5:0.5:3;
ntst = 40;
ncol = 4;
amp = 1e-3;
results = repmat(struct('Q',[], 'gHopf',[], 'l1',[], ...
    'cycleMultipliers',[], 'cycleStable',false), numel(qValues), 1);

for idx = 1:numel(qValues)
    Q = qValues(idx);
    fprintf('\n=== Processing Q = %.2f ===\n', Q);
    p = [0.8; Q];
    ap = 1; % continue in g

    try
        [x0, v0] = init_EP_EP(odefile, [0; 0; 0], p, ap);
    catch ME
        warning('Initialization of equilibrium continuation failed at Q = %.2f: %s', Q, ME.message);
        continue;
    end

    optEQ = contset;
    optEQ = contset(optEQ, 'MaxNumPoints', 200);
    optEQ = contset(optEQ, 'Singularities', 1);
    optEQ = contset(optEQ, 'Eigenvalues', 1);

    [xeq, veq, seq, heq, feq] = cont(@equilibrium, x0, v0, optEQ);

    hopfIdx = find(strcmp({seq.label}, 'H'), 1);
    if isempty(hopfIdx)
        fprintf('  No Hopf point detected on the equilibrium curve for Q = %.2f.\n', Q);
        continue;
    end

    hopfStateIdx = seq(hopfIdx).index;
    xHopf = xeq(1:3, hopfStateIdx);
    gHopf = xeq(end, hopfStateIdx);
    l1 = seq(hopfIdx).data.lyapunov;

    results(idx).Q = Q;
    results(idx).gHopf = gHopf;
    results(idx).l1 = l1;

    hopfType = 'subcritical';
    if l1 < 0
        hopfType = 'supercritical';
    end
    fprintf('  Hopf detected at g = %.6f with l1 = %.4e (%s).\n', gHopf, l1, hopfType);

    figure('Name', sprintf('Equilibria Q=%.2f', Q));
    cpl(xeq, veq, seq, [size(xeq,1) 1]);
    xlabel('g'); ylabel('x'); title(sprintf('Equilibria for Q = %.2f', Q));

    pHopf = [gHopf; Q];
    [xlc0, vlc0] = init_H_LC(odefile, xHopf, pHopf, ap, amp, ntst, ncol);

    optLC = contset;
    optLC = contset(optLC, 'MaxNumPoints', 200);
    optLC = contset(optLC, 'Multipliers', 1);
    optLC = contset(optLC, 'IgnoreSingularity', 1);
    optLC = contset(optLC, 'Adapt', 1);

    [xlc, vlc, slc, hlc, flc] = cont(@limitcycle, xlc0, vlc0, optLC);

    figure('Name', sprintf('Limit cycles Q=%.2f', Q));
    cpl(xlc, vlc, slc, [size(xlc,1) 1]);
    xlabel('g'); ylabel('x'); title(sprintf('Limit cycles for Q = %.2f', Q));

    multipliers = extract_cycle_multipliers(slc, flc, 3);
    results(idx).cycleMultipliers = multipliers;

    if isempty(multipliers)
        fprintf('  Floquet multipliers unavailable for Q = %.2f.\n', Q);
        continue;
    end

    fprintf('  Floquet multipliers:');
    for k = 1:numel(multipliers)
        fprintf(' % .4f%+.4fi', real(multipliers(k)), imag(multipliers(k)));
    end
    fprintf('\n');

    multipliersNoTrivial = remove_trivial_multiplier(multipliers);
    if isempty(multipliersNoTrivial)
        fprintf('  Only the trivial Floquet multiplier was found; stability could not be assessed.\n');
        results(idx).cycleStable = false;
        continue;
    end

    isStable = all(abs(multipliersNoTrivial) < 1 - 1e-3);
    results(idx).cycleStable = isStable;
    if isStable
        fprintf('  Limit cycle is numerically stable (non-trivial multipliers inside unit circle).\n');
    else
        fprintf('  Limit cycle is unstable (a non-trivial multiplier lies on/ outside the unit circle).\n');
    end
end

fprintf('\nSummary over Q values:\n');
for idx = 1:numel(results)
    if isempty(results(idx).Q)
        continue;
    end
    cycleVals = remove_trivial_multiplier(results(idx).cycleMultipliers);
    if isempty(cycleVals)
        multMsg = 'no multipliers';
        stabilityMsg = 'cycle status unknown';
    else
        multMsg = sprintf('max |λ| = %.4f', max(abs(cycleVals)));
        stabilityMsg = ternary(results(idx).cycleStable, 'stable cycle', 'unstable cycle');
    end
    fprintf('  Q = %.2f: g_H = %.6f, l1 = %.4e, %s, %s, %s.\n', ...
        results(idx).Q, results(idx).gHopf, results(idx).l1, ...
        ternary(results(idx).l1 < 0, 'supercritical Hopf', 'subcritical Hopf'), multMsg, stabilityMsg);
end

end

% -------------------------------------------------------------------------
function multipliers = extract_cycle_multipliers(slc, flc, nphase)
multipliers = [];
if ~isempty(slc)
    for idx = numel(slc):-1:1
        if isfield(slc(idx), 'data') && isfield(slc(idx).data, 'multipliers') && ...
                ~isempty(slc(idx).data.multipliers)
            multipliers = slc(idx).data.multipliers(:);
            break;
        end
    end
end

if isempty(multipliers) && ~isempty(flc) && size(flc,1) >= nphase
    multipliers = flc(end-nphase+1:end, end);
end
end

% -------------------------------------------------------------------------
function values = remove_trivial_multiplier(values)
if isempty(values)
    return;
end
[~, idx] = min(abs(values - 1));
values(idx) = [];
end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
