from dataclasses import replace
from typing import Tuple

from .specs import SectionProfileRelationSpec


def independent_profile_relations(section_count: int = 4) -> Tuple[SectionProfileRelationSpec, ...]:
    return tuple(SectionProfileRelationSpec() for _ in range(int(section_count)))


def tie_shape(
    relations: Tuple[SectionProfileRelationSpec, ...],
    target_index: int,
    source_index: int,
) -> Tuple[SectionProfileRelationSpec, ...]:
    updated = list(relations)
    updated[target_index] = replace(updated[target_index], shape_source_index=int(source_index))
    return tuple(updated)


def tie_cst(
    relations: Tuple[SectionProfileRelationSpec, ...],
    target_index: int,
    source_index: int,
) -> Tuple[SectionProfileRelationSpec, ...]:
    updated = list(relations)
    updated[target_index] = replace(
        updated[target_index],
        upper_source_index=int(source_index),
        lower_source_index=int(source_index),
    )
    return tuple(updated)


def tie_upper(
    relations: Tuple[SectionProfileRelationSpec, ...],
    target_index: int,
    source_index: int,
) -> Tuple[SectionProfileRelationSpec, ...]:
    updated = list(relations)
    updated[target_index] = replace(updated[target_index], upper_source_index=int(source_index))
    return tuple(updated)


def tie_lower(
    relations: Tuple[SectionProfileRelationSpec, ...],
    target_index: int,
    source_index: int,
) -> Tuple[SectionProfileRelationSpec, ...]:
    updated = list(relations)
    updated[target_index] = replace(updated[target_index], lower_source_index=int(source_index))
    return tuple(updated)


def tie_te_thickness(
    relations: Tuple[SectionProfileRelationSpec, ...],
    target_index: int,
    source_index: int,
) -> Tuple[SectionProfileRelationSpec, ...]:
    updated = list(relations)
    updated[target_index] = replace(updated[target_index], te_source_index=int(source_index))
    return tuple(updated)
