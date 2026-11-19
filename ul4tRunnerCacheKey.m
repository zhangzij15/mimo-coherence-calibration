function digest = ul4tRunnerCacheKey(cfg,sourceFingerprint)
%UL4TRUNNERCACHEKEY SHA256 of complete, canonically ordered configuration.
payload=struct('schema','matlab5g_coherence_campaign_v3', ...
    'source_sha256',sourceFingerprint,'configuration',canonical(cfg));
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(uint8(unicode2native(jsonencode(payload),'UTF-8')));
raw=typecast(md.digest(),'uint8');
digest=lower(reshape(dec2hex(raw,2).',1,[]));
end

function x=canonical(x)
if isstruct(x)
    x=orderfields(x); fields=fieldnames(x);
    for k=1:numel(x)
        for j=1:numel(fields), x(k).(fields{j})=canonical(x(k).(fields{j})); end
    end
elseif iscell(x)
    for k=1:numel(x), x{k}=canonical(x{k}); end
end
end
