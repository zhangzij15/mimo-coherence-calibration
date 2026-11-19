function [model,scene] = ul4tNativeCDLC(cfg,pattern,schedule,mc,ue)
%UL4TNATIVECDLC Keep the original CDL-C link-adapter entry point available.
assert(strcmpi(cfg.cdl.toolbox_profile,'CDL-C'),'UL256:NativeCDLProfile', ...
    'Use ul4tNativeCDL for a profile other than CDL-C.');
[model,scene] = ul4tNativeCDL(cfg,pattern,schedule,mc,ue);
model.diagnostics.backend = 'base-cdl-c';
end
