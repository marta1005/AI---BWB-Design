from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, Tuple

from .laws import InterpolationSpec, ScalarLawSpec
from .sections import SectionAnchorSpec


class ComponentKind(str, Enum):
    GENERIC_LOFT = "generic_loft"
    LIFTING_SURFACE = "lifting_surface"
    FUSELAGE = "fuselage"
    NACELLE = "nacelle"


@dataclass(frozen=True)
class LoftedComponentSpec:
    component_id: str
    kind: ComponentKind
    sections: Tuple[SectionAnchorSpec, ...]
    section_interpolation: InterpolationSpec = InterpolationSpec()
    spine_interpolation: InterpolationSpec = InterpolationSpec()
    scalar_laws: Tuple[ScalarLawSpec, ...] = ()
    mirrored: bool = False
    metadata: Dict[str, str] = field(default_factory=dict)

    def validate(self) -> None:
        if not self.component_id:
            raise ValueError("component_id must be non-empty")
        if len(self.sections) < 2:
            raise ValueError(
                f"component {self.component_id!r} must define at least 2 anchor sections"
            )
        self.section_interpolation.validate()
        self.spine_interpolation.validate()

        section_ids = set()
        section_etas = []
        for section in self.sections:
            section.validate()
            if section.section_id in section_ids:
                raise ValueError(
                    f"component {self.component_id!r} has duplicated section_id {section.section_id!r}"
                )
            section_ids.add(section.section_id)
            section_etas.append(float(section.eta))

        if any(left >= right for left, right in zip(section_etas[:-1], section_etas[1:])):
            raise ValueError(
                f"component {self.component_id!r} sections must be strictly increasing in eta, "
                f"got {tuple(section_etas)}"
            )

        law_names = set()
        for law in self.scalar_laws:
            law.validate()
            if law.name in law_names:
                raise ValueError(
                    f"component {self.component_id!r} has duplicated scalar law {law.name!r}"
                )
            law_names.add(law.name)

    def section_ids(self) -> Tuple[str, ...]:
        self.validate()
        return tuple(section.section_id for section in self.sections)

    def profile_ids(self) -> Tuple[str, ...]:
        self.validate()
        return tuple(section.profile_id for section in self.sections)

    def law_names(self) -> Tuple[str, ...]:
        self.validate()
        return tuple(law.name for law in self.scalar_laws)

    def scalar_law_map(self) -> Dict[str, ScalarLawSpec]:
        self.validate()
        return {law.name: law for law in self.scalar_laws}
