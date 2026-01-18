<script lang="ts">
	import { onMount } from 'svelte';
	import {
		getMitoMasterStatus,
		getMitoMasterStats,
		getVariantAtPosition,
		searchVariants,
		getGenes,
		type MitoMasterStatus,
		type MitoMasterStats,
		type VariantDetailResponse,
		type VariantFrequency,
		type GeneInfo
	} from '$lib/api';

	// State
	let status: MitoMasterStatus | null = $state(null);
	let stats: MitoMasterStats | null = $state(null);
	let genes: GeneInfo[] = $state([]);
	let loading = $state(true);
	let error = $state('');

	// Search state
	let searchMode = $state<'position' | 'range' | 'gene'>('position');
	let positionInput = $state('');
	let rangeStart = $state('');
	let rangeEnd = $state('');
	let selectedGene = $state('');
	let minFrequency = $state('');

	// Results
	let variantResult: VariantDetailResponse | null = $state(null);
	let searchResults: VariantFrequency[] = $state([]);
	let searchLoading = $state(false);
	let searchError = $state('');

	onMount(async () => {
		try {
			status = await getMitoMasterStatus();

			if (status.status === 'ok') {
				const [statsRes, genesRes] = await Promise.all([getMitoMasterStats(), getGenes()]);
				stats = statsRes;
				genes = genesRes;
			}
		} catch (e) {
			error = e instanceof Error ? e.message : 'Failed to load MitoMaster status';
		} finally {
			loading = false;
		}
	});

	async function handleSearch() {
		searchError = '';
		variantResult = null;
		searchResults = [];
		searchLoading = true;

		try {
			if (searchMode === 'position') {
				const pos = parseInt(positionInput);
				if (isNaN(pos) || pos < 1 || pos > 16569) {
					throw new Error('Position must be between 1 and 16569');
				}
				variantResult = await getVariantAtPosition(pos);
			} else if (searchMode === 'range') {
				const params: any = { limit: 100 };
				if (rangeStart) params.position_start = parseInt(rangeStart);
				if (rangeEnd) params.position_end = parseInt(rangeEnd);
				if (minFrequency) params.min_frequency = parseFloat(minFrequency);

				const res = await searchVariants(params);
				searchResults = res.variants;
			} else if (searchMode === 'gene') {
				if (!selectedGene) {
					throw new Error('Please select a gene');
				}
				const params: any = { gene: selectedGene, limit: 100 };
				if (minFrequency) params.min_frequency = parseFloat(minFrequency);

				const res = await searchVariants(params);
				searchResults = res.variants;
			}
		} catch (e) {
			searchError = e instanceof Error ? e.message : 'Search failed';
		} finally {
			searchLoading = false;
		}
	}

	function formatPercent(value: number): string {
		return (value * 100).toFixed(2) + '%';
	}

	function formatNumber(value: number): string {
		return value.toLocaleString();
	}

	function getGeneTypeColor(type: string): string {
		switch (type) {
			case 'protein_coding':
				return 'bg-blue-100 text-blue-800';
			case 'tRNA':
				return 'bg-green-100 text-green-800';
			case 'rRNA':
				return 'bg-purple-100 text-purple-800';
			case 'control_region':
				return 'bg-amber-100 text-amber-800';
			default:
				return 'bg-slate-100 text-slate-800';
		}
	}
</script>

<svelte:head>
	<title>MitoMaster - Variant Database</title>
</svelte:head>

<main class="max-w-6xl mx-auto px-4 py-8">
	{#if loading}
		<div class="flex items-center justify-center py-12">
			<svg class="w-8 h-8 animate-spin text-blue-600" viewBox="0 0 24 24" fill="none">
				<circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"
				></circle>
				<path
					class="opacity-75"
					fill="currentColor"
					d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
				></path>
			</svg>
		</div>
	{:else if error}
		<div class="bg-red-50 border border-red-200 rounded-lg p-6 text-red-700">
			<h2 class="font-semibold mb-2">Error</h2>
			<p>{error}</p>
		</div>
	{:else if status?.status !== 'ok'}
		<!-- Database not configured -->
		<div class="bg-amber-50 border border-amber-200 rounded-lg p-6">
			<h2 class="text-lg font-semibold text-amber-800 mb-2">MitoMaster Database Not Available</h2>
			<p class="text-amber-700 mb-4">
				{status?.message || 'The MitoMaster database has not been built yet.'}
			</p>
			<div class="bg-white rounded p-4 font-mono text-sm">
				<p class="text-slate-600 mb-2"># Build the database from NCBI:</p>
				<p class="text-slate-800">haplogrep3 mitomaster-build -e your@email.com</p>
			</div>
		</div>
	{:else}
		<!-- Stats Cards -->
		{#if stats}
			<div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
				<div class="bg-white rounded-lg border border-slate-200 p-4">
					<p class="text-xs uppercase tracking-wide text-slate-500">Total Genomes</p>
					<p class="text-2xl font-bold text-slate-800">{formatNumber(stats.genome_count)}</p>
				</div>
				<div class="bg-white rounded-lg border border-slate-200 p-4">
					<p class="text-xs uppercase tracking-wide text-slate-500">Human Genomes</p>
					<p class="text-2xl font-bold text-slate-800">{formatNumber(stats.human_genome_count)}</p>
				</div>
				<div class="bg-white rounded-lg border border-slate-200 p-4">
					<p class="text-xs uppercase tracking-wide text-slate-500">Total Variants</p>
					<p class="text-2xl font-bold text-slate-800">{formatNumber(stats.variant_count)}</p>
				</div>
				<div class="bg-white rounded-lg border border-slate-200 p-4">
					<p class="text-xs uppercase tracking-wide text-slate-500">Unique Variants</p>
					<p class="text-2xl font-bold text-slate-800">{formatNumber(stats.unique_variant_count)}</p>
				</div>
			</div>
		{/if}

		<!-- Search Card -->
		<div class="bg-white rounded-xl shadow-sm border border-slate-200 p-6 mb-6">
			<h2 class="text-lg font-semibold text-slate-800 mb-4 flex items-center gap-2">
				<svg
					class="w-5 h-5 text-blue-600"
					viewBox="0 0 24 24"
					fill="none"
					stroke="currentColor"
					stroke-width="2"
				>
					<circle cx="11" cy="11" r="8" />
					<path d="m21 21-4.35-4.35" />
				</svg>
				Search Variants
			</h2>

			<!-- Search Mode Tabs -->
			<div class="flex gap-2 mb-4">
				<button
					onclick={() => (searchMode = 'position')}
					class="px-4 py-2 rounded-lg text-sm font-medium transition-colors
						{searchMode === 'position'
						? 'bg-blue-600 text-white'
						: 'bg-slate-100 text-slate-700 hover:bg-slate-200'}"
				>
					Position
				</button>
				<button
					onclick={() => (searchMode = 'range')}
					class="px-4 py-2 rounded-lg text-sm font-medium transition-colors
						{searchMode === 'range'
						? 'bg-blue-600 text-white'
						: 'bg-slate-100 text-slate-700 hover:bg-slate-200'}"
				>
					Range
				</button>
				<button
					onclick={() => (searchMode = 'gene')}
					class="px-4 py-2 rounded-lg text-sm font-medium transition-colors
						{searchMode === 'gene'
						? 'bg-blue-600 text-white'
						: 'bg-slate-100 text-slate-700 hover:bg-slate-200'}"
				>
					Gene
				</button>
			</div>

			<!-- Search Inputs -->
			<div class="space-y-4">
				{#if searchMode === 'position'}
					<div>
						<label for="position" class="block text-sm font-medium text-slate-700 mb-1"
							>mtDNA Position (1-16569)</label
						>
						<input
							type="number"
							id="position"
							bind:value={positionInput}
							min="1"
							max="16569"
							placeholder="e.g., 3243"
							class="block w-full h-10 px-3 bg-white border border-slate-300 rounded-md text-sm shadow-sm
								focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
						/>
					</div>
				{:else if searchMode === 'range'}
					<div class="grid grid-cols-2 gap-4">
						<div>
							<label for="range-start" class="block text-sm font-medium text-slate-700 mb-1"
								>Start Position</label
							>
							<input
								type="number"
								id="range-start"
								bind:value={rangeStart}
								min="1"
								max="16569"
								placeholder="1"
								class="block w-full h-10 px-3 bg-white border border-slate-300 rounded-md text-sm shadow-sm
									focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
							/>
						</div>
						<div>
							<label for="range-end" class="block text-sm font-medium text-slate-700 mb-1"
								>End Position</label
							>
							<input
								type="number"
								id="range-end"
								bind:value={rangeEnd}
								min="1"
								max="16569"
								placeholder="16569"
								class="block w-full h-10 px-3 bg-white border border-slate-300 rounded-md text-sm shadow-sm
									focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
							/>
						</div>
					</div>
				{:else if searchMode === 'gene'}
					<div>
						<label for="gene" class="block text-sm font-medium text-slate-700 mb-1">Gene</label>
						<select
							id="gene"
							bind:value={selectedGene}
							class="block w-full h-10 px-3 bg-white border border-slate-300 rounded-md text-sm shadow-sm
								focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
						>
							<option value="">Select a gene...</option>
							{#each genes as gene (gene.name)}
								<option value={gene.name}
									>{gene.name} ({gene.start}-{gene.end}) - {gene.product}</option
								>
							{/each}
						</select>
					</div>
				{/if}

				{#if searchMode !== 'position'}
					<div>
						<label for="min-freq" class="block text-sm font-medium text-slate-700 mb-1"
							>Min Frequency (optional)</label
						>
						<input
							type="number"
							id="min-freq"
							bind:value={minFrequency}
							min="0"
							max="1"
							step="0.01"
							placeholder="e.g., 0.01"
							class="block w-full h-10 px-3 bg-white border border-slate-300 rounded-md text-sm shadow-sm
								focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
						/>
					</div>
				{/if}

				<button
					onclick={handleSearch}
					disabled={searchLoading}
					class="w-full py-3 px-4 bg-blue-600 hover:bg-blue-700 disabled:bg-slate-300 disabled:cursor-not-allowed text-white font-medium rounded-lg transition-colors flex items-center justify-center gap-2"
				>
					{#if searchLoading}
						<svg class="w-5 h-5 animate-spin" viewBox="0 0 24 24" fill="none">
							<circle
								class="opacity-25"
								cx="12"
								cy="12"
								r="10"
								stroke="currentColor"
								stroke-width="4"
							></circle>
							<path
								class="opacity-75"
								fill="currentColor"
								d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
							></path>
						</svg>
						Searching...
					{:else}
						Search
					{/if}
				</button>

				{#if searchError}
					<div class="p-4 bg-red-50 border border-red-200 rounded-lg text-red-700 text-sm">
						{searchError}
					</div>
				{/if}
			</div>
		</div>

		<!-- Position Search Result -->
		{#if variantResult}
			<div class="bg-white rounded-xl shadow-sm border border-slate-200 p-6 mb-6">
				<h2 class="text-lg font-semibold text-slate-800 mb-4">
					Position {variantResult.variant.position}
					{#if variantResult.variant.locus}
						<span class="text-slate-500 font-normal">({variantResult.variant.locus})</span>
					{/if}
				</h2>

				<div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
					<div class="bg-slate-50 rounded-lg p-4">
						<p class="text-xs uppercase tracking-wide text-slate-500">Variant</p>
						<p class="text-xl font-mono font-bold text-slate-800">
							{variantResult.variant.ref}>{variantResult.variant.alt}
						</p>
					</div>
					<div class="bg-slate-50 rounded-lg p-4">
						<p class="text-xs uppercase tracking-wide text-slate-500">Total Count</p>
						<p class="text-xl font-bold text-slate-800">
							{formatNumber(variantResult.variant.total_count)}
						</p>
					</div>
					<div class="bg-slate-50 rounded-lg p-4">
						<p class="text-xs uppercase tracking-wide text-slate-500">Total Frequency</p>
						<p class="text-xl font-bold text-slate-800">
							{formatPercent(variantResult.variant.total_frequency)}
						</p>
					</div>
					<div class="bg-slate-50 rounded-lg p-4">
						<p class="text-xs uppercase tracking-wide text-slate-500">Human Frequency</p>
						<p class="text-xl font-bold text-slate-800">
							{formatPercent(variantResult.variant.human_frequency)}
						</p>
					</div>
				</div>

				{#if variantResult.variant.gene_type}
					<div class="mb-4">
						<span class="inline-block px-2 py-1 rounded text-xs font-medium {getGeneTypeColor(variantResult.variant.gene_type)}">
							{variantResult.variant.gene_type}
						</span>
						{#if variantResult.variant.amino_acid_change}
							<span class="ml-2 font-mono text-sm">{variantResult.variant.amino_acid_change}</span>
						{/if}
					</div>
				{/if}

				{#if variantResult.haplogroup_frequencies.length > 0}
					<h3 class="font-semibold text-slate-700 mb-3">Haplogroup Distribution</h3>
					<div class="overflow-x-auto">
						<table class="w-full text-sm">
							<thead>
								<tr class="border-b border-slate-200">
									<th class="text-left py-2 px-4 font-semibold text-slate-600">Haplogroup</th>
									<th class="text-right py-2 px-4 font-semibold text-slate-600">Count</th>
									<th class="text-right py-2 px-4 font-semibold text-slate-600">Frequency</th>
									<th class="py-2 px-4"></th>
								</tr>
							</thead>
							<tbody>
								{#each variantResult.haplogroup_frequencies as hf (hf.haplogroup)}
									<tr class="border-b border-slate-100">
										<td class="py-2 px-4 font-semibold text-blue-600">{hf.haplogroup}</td>
										<td class="py-2 px-4 text-right">{formatNumber(hf.count)}</td>
										<td class="py-2 px-4 text-right">{formatPercent(hf.frequency)}</td>
										<td class="py-2 px-4">
											<div class="w-24 bg-slate-200 rounded-full h-2">
												<div
													class="bg-blue-600 h-2 rounded-full"
													style="width: {Math.min(hf.frequency * 100, 100)}%"
												></div>
											</div>
										</td>
									</tr>
								{/each}
							</tbody>
						</table>
					</div>
				{/if}
			</div>
		{/if}

		<!-- Range/Gene Search Results -->
		{#if searchResults.length > 0}
			<div class="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
				<h2 class="text-lg font-semibold text-slate-800 mb-4">
					Search Results ({searchResults.length} variants)
				</h2>

				<div class="overflow-x-auto">
					<table class="w-full text-sm">
						<thead>
							<tr class="border-b border-slate-200">
								<th class="text-left py-2 px-4 font-semibold text-slate-600">Position</th>
								<th class="text-left py-2 px-4 font-semibold text-slate-600">Variant</th>
								<th class="text-left py-2 px-4 font-semibold text-slate-600">Gene</th>
								<th class="text-right py-2 px-4 font-semibold text-slate-600">Count</th>
								<th class="text-right py-2 px-4 font-semibold text-slate-600">Frequency</th>
								<th class="text-left py-2 px-4 font-semibold text-slate-600">AA Change</th>
							</tr>
						</thead>
						<tbody>
							{#each searchResults as v (v.position + v.ref + v.alt)}
								<tr class="border-b border-slate-100 hover:bg-slate-50">
									<td class="py-2 px-4">
										<button
											onclick={() => {
												positionInput = v.position.toString();
												searchMode = 'position';
												handleSearch();
											}}
											class="text-blue-600 hover:underline"
										>
											{v.position}
										</button>
									</td>
									<td class="py-2 px-4 font-mono">{v.ref}>{v.alt}</td>
									<td class="py-2 px-4">
										{#if v.locus}
											<span class="inline-block px-2 py-0.5 rounded text-xs font-medium {getGeneTypeColor(v.gene_type || '')}">
												{v.locus}
											</span>
										{/if}
									</td>
									<td class="py-2 px-4 text-right">{formatNumber(v.total_count)}</td>
									<td class="py-2 px-4 text-right">{formatPercent(v.total_frequency)}</td>
									<td class="py-2 px-4 font-mono text-xs">{v.amino_acid_change || ''}</td>
								</tr>
							{/each}
						</tbody>
					</table>
				</div>
			</div>
		{/if}

		<!-- Top Haplogroups -->
		{#if stats && stats.top_haplogroups.length > 0 && !variantResult && searchResults.length === 0}
			<div class="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
				<h2 class="text-lg font-semibold text-slate-800 mb-4">Top Haplogroups in Database</h2>
				<div class="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-3">
					{#each stats.top_haplogroups as hg (hg.haplogroup)}
						<div class="bg-slate-50 rounded-lg p-3 text-center">
							<p class="font-semibold text-blue-600">{hg.haplogroup}</p>
							<p class="text-sm text-slate-500">{formatNumber(hg.count)}</p>
						</div>
					{/each}
				</div>
			</div>
		{/if}
	{/if}
</main>
