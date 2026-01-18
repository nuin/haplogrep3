// API types matching FastAPI response models

export interface TreeInfo {
	id: string;
	name: string;
	version: string;
}

export interface TreesResponse {
	trees: TreeInfo[];
}

export interface DistanceInfo {
	id: string;
	name: string;
	default: boolean;
}

export interface DistancesResponse {
	distances: DistanceInfo[];
}

export interface ClassificationResultResponse {
	haplogroup: string;
	quality: number;
	found_count: number;
	missing_count: number;
	remaining_count: number;
	found_polymorphisms: string[];
	missing_polymorphisms: string[];
	remaining_polymorphisms: string[];
}

export interface SampleResultResponse {
	id: string;
	haplogroup: string;
	quality: number;
	range: string;
	polymorphisms: string[];
	top_result: ClassificationResultResponse;
	other_results: ClassificationResultResponse[];
}

export interface ClassifyResponse {
	job_id: string;
	tree: string;
	distance: string;
	samples_count: number;
	results: SampleResultResponse[];
}

export interface HealthResponse {
	status: string;
	version: string;
}

// API client functions

const API_BASE = '/api';

export async function checkHealth(): Promise<HealthResponse> {
	const response = await fetch(`${API_BASE}/health`);
	if (!response.ok) {
		throw new Error(`Health check failed: ${response.statusText}`);
	}
	return response.json();
}

export async function getTrees(): Promise<TreesResponse> {
	const response = await fetch(`${API_BASE}/trees`);
	if (!response.ok) {
		throw new Error(`Failed to fetch trees: ${response.statusText}`);
	}
	return response.json();
}

export async function getDistances(): Promise<DistancesResponse> {
	const response = await fetch(`${API_BASE}/distances`);
	if (!response.ok) {
		throw new Error(`Failed to fetch distances: ${response.statusText}`);
	}
	return response.json();
}

export async function classify(
	file: File,
	tree: string,
	distance: string = 'kulczynski',
	hits: number = 1,
	hetLevel: number = 0.9,
	chip: boolean = false
): Promise<ClassifyResponse> {
	const formData = new FormData();
	formData.append('file', file);
	formData.append('tree', tree);
	formData.append('distance', distance);
	formData.append('hits', hits.toString());
	formData.append('het_level', hetLevel.toString());
	formData.append('chip', chip.toString());

	const response = await fetch(`${API_BASE}/classify`, {
		method: 'POST',
		body: formData
	});

	if (!response.ok) {
		const error = await response.json();
		throw new Error(error.detail || `Classification failed: ${response.statusText}`);
	}

	return response.json();
}

// =============================================================================
// MitoMaster Types
// =============================================================================

export interface VariantFrequency {
	position: number;
	ref: string;
	alt: string;
	total_count: number;
	total_frequency: number;
	human_count: number;
	human_frequency: number;
	locus: string | null;
	gene_type: string | null;
	amino_acid_change: string | null;
}

export interface HaplogroupFrequency {
	haplogroup: string;
	count: number;
	frequency: number;
}

export interface VariantDetailResponse {
	variant: VariantFrequency;
	haplogroup_frequencies: HaplogroupFrequency[];
}

export interface GeneInfo {
	name: string;
	start: number;
	end: number;
	type: string;
	strand: string;
	product: string;
}

export interface MitoMasterStats {
	genome_count: number;
	human_genome_count: number;
	variant_count: number;
	unique_variant_count: number;
	top_haplogroups: { haplogroup: string; count: number }[];
}

export interface MitoMasterStatus {
	status: string;
	message?: string;
	database_path: string | null;
	genome_count?: number;
	variant_count?: number;
}

export interface SearchVariantsResponse {
	count: number;
	variants: VariantFrequency[];
}

// =============================================================================
// MitoMaster API Functions
// =============================================================================

export async function getMitoMasterStatus(): Promise<MitoMasterStatus> {
	const response = await fetch(`${API_BASE}/mitomaster/status`);
	if (!response.ok) {
		throw new Error(`Failed to get MitoMaster status: ${response.statusText}`);
	}
	return response.json();
}

export async function getMitoMasterStats(): Promise<MitoMasterStats> {
	const response = await fetch(`${API_BASE}/mitomaster/stats`);
	if (!response.ok) {
		throw new Error(`Failed to get MitoMaster stats: ${response.statusText}`);
	}
	return response.json();
}

export async function getVariantAtPosition(
	position: number,
	ref?: string,
	alt?: string
): Promise<VariantDetailResponse> {
	const params = new URLSearchParams();
	if (ref) params.append('ref', ref);
	if (alt) params.append('alt', alt);
	const queryString = params.toString() ? `?${params.toString()}` : '';

	const response = await fetch(`${API_BASE}/mitomaster/variant/${position}${queryString}`);
	if (!response.ok) {
		if (response.status === 404) {
			throw new Error(`No variants found at position ${position}`);
		}
		throw new Error(`Failed to get variant: ${response.statusText}`);
	}
	return response.json();
}

export async function searchVariants(params: {
	position_start?: number;
	position_end?: number;
	gene?: string;
	min_frequency?: number;
	limit?: number;
}): Promise<SearchVariantsResponse> {
	const searchParams = new URLSearchParams();
	if (params.position_start) searchParams.append('position_start', params.position_start.toString());
	if (params.position_end) searchParams.append('position_end', params.position_end.toString());
	if (params.gene) searchParams.append('gene', params.gene);
	if (params.min_frequency) searchParams.append('min_frequency', params.min_frequency.toString());
	if (params.limit) searchParams.append('limit', params.limit.toString());

	const queryString = searchParams.toString() ? `?${searchParams.toString()}` : '';
	const response = await fetch(`${API_BASE}/mitomaster/search${queryString}`);
	if (!response.ok) {
		throw new Error(`Failed to search variants: ${response.statusText}`);
	}
	return response.json();
}

export async function getGenes(): Promise<GeneInfo[]> {
	const response = await fetch(`${API_BASE}/mitomaster/genes`);
	if (!response.ok) {
		throw new Error(`Failed to get genes: ${response.statusText}`);
	}
	return response.json();
}
