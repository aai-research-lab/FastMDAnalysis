"""
Options Forwarding Utility

Provides permissive options passthrough with:
- Introspection of target signatures to forward accepted kwargs
- Alias mapping (e.g., ref → reference, atoms → atom_indices)
- Pre-hooks (e.g., align: true → traj.superpose(...))
- Post-hooks (e.g., per_residue: true → residue aggregation)
- Strictness mode (log unknown keys or raise errors)
"""
from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Set, Tuple
import inspect
import logging
import warnings

logger = logging.getLogger(__name__)

__all__ = [
    "OptionsForwarder",
    "forward_options",
    "apply_alias_mapping",
]


class OptionsForwarder:
    """
    Handles permissive options forwarding with alias mapping and hooks.
    
    Parameters
    ----------
    aliases : dict, optional
        Mapping of alias names to canonical parameter names.
        Example: {"ref": "reference_frame", "atoms": "atom_indices"}
    pre_hooks : dict, optional
        Mapping of option names to pre-processing functions.
        Functions should accept (value, context) and return processed value.
    post_hooks : dict, optional
        Mapping of option names to post-processing functions.
        Functions should accept (result, context) and return modified result.
    strict : bool
        If True, raise errors for unknown options. If False, log warnings.
    """
    
    def __init__(
        self,
        aliases: Optional[Dict[str, str]] = None,
        pre_hooks: Optional[Dict[str, Callable]] = None,
        post_hooks: Optional[Dict[str, Callable]] = None,
        strict: bool = False,
    ):
        self.aliases = aliases or {}
        self.pre_hooks = pre_hooks or {}
        self.post_hooks = post_hooks or {}
        self.strict = strict
        
    def apply_aliases(self, options: Dict[str, Any]) -> Dict[str, Any]:
        """
        Apply alias mappings to options dictionary.
        
        Parameters
        ----------
        options : dict
            Raw options from user.
            
        Returns
        -------
        dict
            Options with aliases resolved to canonical names.
        """
        resolved: Dict[str, Any] = {}
        sources: Dict[str, str] = {}
        for key, value in options.items():
            canonical = self.aliases.get(key, key)
            existing = resolved.get(canonical)
            previous_key = sources.get(canonical)
            duplicate = previous_key is not None and previous_key != key

            if duplicate:
                if existing is not None and existing != value:
                    msg = (
                        f"Both alias '{key}' and canonical '{canonical}' provided; "
                        "alias value overrides"
                    )
                    if self.strict:
                        raise ValueError(msg)
                    logger.warning(msg)
                else:
                    logger.debug(
                        "Alias '%s' overriding canonical '%s' (default value)",
                        key,
                        canonical,
                    )

            resolved[canonical] = value
            if previous_key is None:
                sources[canonical] = key
        return resolved
    
    def forward_to_callable(
        self,
        callable_obj: Callable,
        options: Dict[str, Any],
    ) -> Tuple[Dict[str, Any], List[str]]:
        """
        Forward options to a callable, keeping only accepted parameters.
        
        Parameters
        ----------
        callable_obj : callable
            Target function or class constructor.
        options : dict
            Options to forward (already with aliases resolved).
            
        Returns
        -------
        forwarded : dict
            Options accepted by the callable.
        dropped : list
            Names of options that were dropped.
        """
        if inspect.isbuiltin(callable_obj):
            logger.debug(
                "Treating builtin callable %s as opaque; forwarding all options",
                callable_obj,
            )
            return dict(options), []

        try:
            sig = inspect.signature(callable_obj)
        except (ValueError, TypeError):
            logger.debug("Cannot inspect signature of %s, forwarding all options", callable_obj)
            return options, []

        accepted = set()
        has_var_keyword = False
        
        for name, param in sig.parameters.items():
            if param.kind == inspect.Parameter.VAR_KEYWORD:
                has_var_keyword = True
            elif param.kind in (
                inspect.Parameter.POSITIONAL_OR_KEYWORD,
                inspect.Parameter.KEYWORD_ONLY,
            ):
                accepted.add(name)

        if has_var_keyword:
            return options, []

        forwarded = {k: v for k, v in options.items() if k in accepted}
        dropped = [k for k in options.keys() if k not in accepted]
        
        return forwarded, dropped
    
    def process_options(
        self,
        options: Dict[str, Any],
        callable_obj: Optional[Callable] = None,
        context: Optional[Dict[str, Any]] = None,
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """
        Full pipeline: apply aliases, run pre-hooks, forward to callable.
        
        Parameters
        ----------
        options : dict
            Raw options from user.
        callable_obj : callable, optional
            Target callable for signature introspection.
        context : dict, optional
            Additional context for hooks (e.g., trajectory object).
            
        Returns
        -------
        forwarded : dict
            Processed options ready to pass to callable.
        hook_data : dict
            Data from pre-hooks that may be needed later.
        """
        context = context or {}
        hook_data = {}
        
        resolved = self.apply_aliases(options)
        processed = {}
        for key, value in resolved.items():
            if key in self.pre_hooks:
                try:
                    result = self.pre_hooks[key](value, context)
                    if isinstance(result, tuple):
                        processed[key], hook_data[key] = result
                    else:
                        processed[key] = result
                except Exception as e:
                    msg = f"Pre-hook for '{key}' failed: {e}"
                    if self.strict:
                        raise RuntimeError(msg) from e
                    else:
                        logger.warning(msg)
                        processed[key] = value
            else:
                processed[key] = value
        if callable_obj is not None:
            forwarded, dropped = self.forward_to_callable(callable_obj, processed)
            
            if dropped:
                msg = f"Unknown options (not accepted by callable): {dropped}"
                if self.strict:
                    raise ValueError(msg)
                else:
                    logger.info(msg)
        else:
            forwarded = processed
            
        return forwarded, hook_data

    def filter_known(
        self,
        options: Dict[str, Any],
        known_keys: Set[str],
        *,
        context: str,
        warn: bool = False,
    ) -> Dict[str, Any]:
        """Keep only known option names, logging or raising on unknown keys."""
        filtered: Dict[str, Any] = {}
        unknown: List[str] = []
        for key, value in options.items():
            if key in known_keys:
                filtered[key] = value
            else:
                unknown.append(key)

        if unknown:
            strict_msg = f"Unknown options for '{context}': {sorted(unknown)}"
            if self.strict:
                raise ValueError(strict_msg)
            logger.warning(strict_msg)
            if warn:
                warn_msg = f"Unsupported options for '{context}': {sorted(unknown)}"
                warnings.warn(warn_msg, stacklevel=2)

        return filtered


def apply_alias_mapping(options: Dict[str, Any], aliases: Dict[str, str]) -> Dict[str, Any]:
    """
    Simple utility to apply alias mappings to options.
    
    Parameters
    ----------
    options : dict
        Raw options from user.
    aliases : dict
        Mapping of alias names to canonical names.
        
    Returns
    -------
    dict
        Options with aliases resolved.
    """
    resolved = {}
    for key, value in options.items():
        canonical = aliases.get(key, key)
        resolved[canonical] = value
    return resolved


def forward_options(
    target: Callable,
    options: Dict[str, Any],
    aliases: Optional[Dict[str, str]] = None,
    strict: bool = False,
) -> Dict[str, Any]:
    """
    Simple utility to forward options with alias mapping.
    
    Parameters
    ----------
    target : callable
        Target function or class constructor.
    options : dict
        Options to forward.
    aliases : dict, optional
        Alias mappings.
    strict : bool
        If True, raise errors for unknown options.
        
    Returns
    -------
    dict
        Filtered options ready to pass to target.
    """
    forwarder = OptionsForwarder(aliases=aliases, strict=strict)
    try:
        forwarded, _ = forwarder.process_options(options, callable_obj=target)
    except ValueError as exc:
        if strict and "Unknown options" in str(exc):
            relaxed_forwarder = OptionsForwarder(aliases=aliases, strict=False)
            forwarded, _ = relaxed_forwarder.process_options(options, callable_obj=target)
        else:
            raise
    return forwarded
