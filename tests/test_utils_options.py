import inspect
import pytest

from fastmdanalysis.utils.options import (
    OptionsForwarder,
    apply_alias_mapping,
    forward_options,
)


def test_apply_aliases_alias_overrides_canonical_when_not_strict():
    forwarder = OptionsForwarder(aliases={"selection": "atoms"})
    resolved = forwarder.apply_aliases({"atoms": "CA", "selection": "protein"})
    assert resolved["atoms"] == "protein"


def test_apply_aliases_duplicate_raises_when_strict():
    forwarder = OptionsForwarder(aliases={"selection": "atoms"}, strict=True)
    with pytest.raises(ValueError):
        forwarder.apply_aliases({"atoms": "CA", "selection": "protein"})


def test_forward_to_callable_filters_unknown_kwargs():
    def target(foo, bar=1):
        return foo + bar

    forwarder = OptionsForwarder()
    forwarded, dropped = forwarder.forward_to_callable(target, {"foo": 2, "extra": 5})
    assert forwarded == {"foo": 2}
    assert dropped == ["extra"]


def test_process_options_runs_hooks_and_collects_hook_data():
    def hook(value, context):
        return value * context["factor"], {"original": value}

    forwarder = OptionsForwarder(
        aliases={"alias_scale": "scale"},
        pre_hooks={"scale": hook},
    )

    forwarded, hook_data = forwarder.process_options(
        {"alias_scale": 2, "other": 7, "unused": True},
        callable_obj=lambda scale, other: (scale, other),
        context={"factor": 3},
    )

    assert forwarded == {"scale": 6, "other": 7}
    assert hook_data["scale"] == {"original": 2}


def test_process_options_pre_hook_failure_respects_strict_mode():
    def failing_hook(value, _context):
        raise RuntimeError("boom")

    forwarder_relaxed = OptionsForwarder(pre_hooks={"flag": failing_hook}, strict=False)
    forwarded_relaxed, _ = forwarder_relaxed.process_options({"flag": True})
    assert forwarded_relaxed["flag"] is True

    forwarder_strict = OptionsForwarder(pre_hooks={"flag": failing_hook}, strict=True)
    with pytest.raises(RuntimeError):
        forwarder_strict.process_options({"flag": True})


def test_filter_known_warns_when_requested():
    forwarder = OptionsForwarder(strict=False)
    with pytest.warns(UserWarning):
        filtered = forwarder.filter_known(
            {"keep": 1, "drop": 2},
            {"keep"},
            context="demo",
            warn=True,
        )
    assert filtered == {"keep": 1}


def test_filter_known_strict_raises_on_unknown():
    forwarder = OptionsForwarder(strict=True)
    with pytest.raises(ValueError):
        forwarder.filter_known({"foo": 1}, {"bar"}, context="demo")


def test_apply_alias_mapping_helper_resolves_keys():
    result = apply_alias_mapping({"ref": 1, "atoms": [0, 1]}, {"ref": "reference_frame"})
    assert result == {"reference_frame": 1, "atoms": [0, 1]}


def test_forward_options_falls_back_for_unknown_when_strict():
    def target(foo):
        return foo

    resolved = forward_options(target, {"foo": 1, "unexpected": 7}, strict=True)
    assert resolved == {"foo": 1}


def test_forward_to_callable_handles_builtin_passthrough():
    forwarder = OptionsForwarder()
    forwarded, dropped = forwarder.forward_to_callable(len, {"foo": 1})
    assert forwarded == {"foo": 1}
    assert dropped == []


def test_forward_to_callable_handles_uninspectable_callable(monkeypatch):
    class Callable:
        def __call__(self, *args, **kwargs):
            return args, kwargs

    opaque = Callable()
    def fake_signature(_):
        raise TypeError("cannot inspect")

    monkeypatch.setattr(inspect, "signature", fake_signature)
    forwarder = OptionsForwarder()
    forwarded, dropped = forwarder.forward_to_callable(opaque, {"foo": 1})

    assert forwarded == {"foo": 1}
    assert dropped == []
