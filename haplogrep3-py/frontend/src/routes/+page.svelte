<script lang="ts">
	import {
		getTrees,
		getDistances,
		classify,
		type TreeInfo,
		type DistanceInfo,
		type ClassifyResponse
	} from '$lib/api';
	import { onMount } from 'svelte';
	import { SvelteSet } from 'svelte/reactivity';

	let trees: TreeInfo[] = $state([]);
	let distances: DistanceInfo[] = $state([]);
	let selectedTree = $state('');
	let selectedDistance = $state('kulczynski');
	let hits = $state(1);
	let hetLevel = $state(0.9);
	let chip = $state(false);

	let file: File | null = $state(null);
	let loading = $state(false);
	let error = $state('');
	let results: ClassifyResponse | null = $state(null);
	let expandedSamples = new SvelteSet<string>();
	let dragOver = $state(false);

	onMount(async () => {
		try {
			const [treesRes, distancesRes] = await Promise.all([getTrees(), getDistances()]);
			trees = treesRes.trees;
			distances = distancesRes.distances;

			if (trees.length > 0) {
				// Prefer phylotree-rcrs v17.2 as default
				const preferredTree = trees.find((t) => t.id === 'phylotree-rcrs' && t.version === '17.2');
				selectedTree = preferredTree ? preferredTree.id : trees[0].id;
			}
			const defaultDistance = distances.find((d) => d.default);
			if (defaultDistance) {
				selectedDistance = defaultDistance.id;
			}
		} catch (e) {
			error = e instanceof Error ? e.message : 'Failed to load configuration';
		}
	});

	function handleFileChange(event: Event) {
		const target = event.target as HTMLInputElement;
		if (target.files && target.files.length > 0) {
			file = target.files[0];
			results = null;
			error = '';
		}
	}

	function handleDrop(event: DragEvent) {
		event.preventDefault();
		dragOver = false;
		if (event.dataTransfer?.files && event.dataTransfer.files.length > 0) {
			file = event.dataTransfer.files[0];
			results = null;
			error = '';
		}
	}

	function handleDragOver(event: DragEvent) {
		event.preventDefault();
		dragOver = true;
	}

	function handleDragLeave() {
		dragOver = false;
	}

	async function handleSubmit() {
		if (!file || !selectedTree) {
			error = 'Please select a file and tree';
			return;
		}

		loading = true;
		error = '';
		results = null;

		try {
			results = await classify(file, selectedTree, selectedDistance, hits, hetLevel, chip);
		} catch (e) {
			error = e instanceof Error ? e.message : 'Classification failed';
		} finally {
			loading = false;
		}
	}

	function toggleExpanded(sampleId: string) {
		if (expandedSamples.has(sampleId)) {
			expandedSamples.delete(sampleId);
		} else {
			expandedSamples.add(sampleId);
		}
	}

	function getQualityClasses(quality: number): string {
		if (quality >= 0.9) return 'bg-emerald-500';
		if (quality >= 0.7) return 'bg-amber-500';
		return 'bg-red-500';
	}

	function formatPercent(value: number): string {
		return (value * 100).toFixed(1) + '%';
	}
</script>

<svelte:head>
	<title>Haplogrep3 - mtDNA Classification</title>
</svelte:head>

<!-- Main Content -->
<main class="max-w-6xl mx-auto px-4 py-8">
		<!-- Upload Card -->
		<div class="bg-white rounded-xl shadow-sm border border-slate-200 p-6 mb-6">
			<h2 class="text-lg font-semibold text-slate-800 mb-6 flex items-center gap-2">
				<svg class="w-5 h-5 text-blue-600" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
					<path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4" />
					<polyline points="17,8 12,3 7,8" />
					<line x1="12" y1="3" x2="12" y2="15" />
				</svg>
				Upload Sample File
			</h2>

			<!-- Dropzone -->
			<div
				class="relative border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition-all mb-6
					{dragOver ? 'border-blue-500 bg-blue-50' : file ? 'border-emerald-500 bg-emerald-50' : 'border-slate-300 hover:border-blue-400 hover:bg-slate-50'}"
				ondrop={handleDrop}
				ondragover={handleDragOver}
				ondragleave={handleDragLeave}
				role="button"
				tabindex="0"
			>
				<input
					type="file"
					class="absolute inset-0 w-full h-full opacity-0 cursor-pointer"
					accept=".vcf,.vcf.gz,.fasta,.fa,.fasta.gz,.fa.gz,.tsv,.txt,.csv"
					onchange={handleFileChange}
				/>
				{#if file}
					<div class="flex items-center justify-center gap-2 text-emerald-600">
						<svg class="w-6 h-6" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
							<path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z" />
							<polyline points="14,2 14,8 20,8" />
						</svg>
						<span class="font-medium">{file.name}</span>
					</div>
				{:else}
					<svg class="w-12 h-12 mx-auto text-slate-400 mb-3" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
						<path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4" />
						<polyline points="17,8 12,3 7,8" />
						<line x1="12" y1="3" x2="12" y2="15" />
					</svg>
					<p class="font-medium text-slate-700">Drop VCF, FASTA, or TSV file here</p>
					<p class="text-sm text-slate-500 mt-1">or click to browse</p>
				{/if}
			</div>

			<!-- Options Grid -->
			<div class="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4 mb-6">
				<div class="space-y-1.5">
					<label for="tree-select" class="block text-sm font-medium text-slate-700">Phylogenetic Tree</label>
					<select
						id="tree-select"
						bind:value={selectedTree}
						disabled={trees.length === 0}
						class="block w-full h-10 px-3 bg-white border border-slate-300 rounded-md text-sm shadow-sm
							focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500
							disabled:bg-slate-100 disabled:cursor-not-allowed"
					>
						{#if trees.length === 0}
							<option value="">No trees installed</option>
						{:else}
							{#each trees as tree (tree.id)}
								<option value={tree.id}>{tree.name} v{tree.version}</option>
							{/each}
						{/if}
					</select>
				</div>

				<div class="space-y-1.5">
					<label for="distance-select" class="block text-sm font-medium text-slate-700">Distance Metric</label>
					<select
						id="distance-select"
						bind:value={selectedDistance}
						class="block w-full h-10 px-3 bg-white border border-slate-300 rounded-md text-sm shadow-sm
							focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
					>
						{#each distances as distance (distance.id)}
							<option value={distance.id}>{distance.name}</option>
						{/each}
					</select>
				</div>

				<div class="space-y-1.5">
					<label for="hits-input" class="block text-sm font-medium text-slate-700">Top Hits</label>
					<input
						type="number"
						id="hits-input"
						bind:value={hits}
						min="1"
						max="10"
						class="block w-full h-10 px-3 bg-white border border-slate-300 rounded-md text-sm shadow-sm
							focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
					/>
				</div>

				<div class="space-y-1.5">
					<label for="het-level-input" class="block text-sm font-medium text-slate-700">Heteroplasmy Level</label>
					<input
						type="number"
						id="het-level-input"
						bind:value={hetLevel}
						min="0"
						max="1"
						step="0.05"
						class="block w-full h-10 px-3 bg-white border border-slate-300 rounded-md text-sm shadow-sm
							focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
					/>
				</div>
			</div>

			<!-- Checkbox -->
			<div class="flex items-center mb-6">
				<input
					type="checkbox"
					id="chip-checkbox"
					bind:checked={chip}
					class="h-4 w-4 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
				/>
				<label for="chip-checkbox" class="ml-2 text-sm text-slate-700 cursor-pointer">
					Genotyping chip data
				</label>
			</div>

			<!-- Submit Button -->
			<button
				onclick={handleSubmit}
				disabled={loading || !file || !selectedTree}
				class="w-full py-3 px-4 bg-blue-600 hover:bg-blue-700 disabled:bg-slate-300 disabled:cursor-not-allowed text-white font-medium rounded-lg transition-colors flex items-center justify-center gap-2"
			>
				{#if loading}
					<svg class="w-5 h-5 animate-spin" viewBox="0 0 24 24" fill="none">
						<circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle>
						<path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
					</svg>
					Classifying...
				{:else}
					<svg class="w-5 h-5" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
						<polyline points="16,3 21,3 21,8" />
						<line x1="4" y1="20" x2="21" y2="3" />
						<polyline points="21,16 21,21 16,21" />
						<line x1="15" y1="15" x2="21" y2="21" />
						<line x1="4" y1="4" x2="9" y2="9" />
					</svg>
					Classify Samples
				{/if}
			</button>

			{#if error}
				<div class="mt-4 p-4 bg-red-50 border border-red-200 rounded-lg text-red-700 text-sm">
					{error}
				</div>
			{/if}
		</div>

		<!-- Results Card -->
		{#if results}
			<div class="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
				<h2 class="text-lg font-semibold text-slate-800 mb-6 flex items-center gap-2">
					<svg class="w-5 h-5 text-blue-600" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
						<path d="M22 12h-4l-3 9L9 3l-3 9H2" />
					</svg>
					Classification Results
				</h2>

				<!-- Meta Info -->
				<div class="flex flex-wrap gap-6 p-4 bg-slate-50 rounded-lg mb-6">
					<div>
						<p class="text-xs uppercase tracking-wide text-slate-500">Samples</p>
						<p class="font-semibold text-slate-800">{results.samples_count}</p>
					</div>
					<div>
						<p class="text-xs uppercase tracking-wide text-slate-500">Tree</p>
						<p class="font-semibold text-slate-800">{results.tree}</p>
					</div>
					<div>
						<p class="text-xs uppercase tracking-wide text-slate-500">Distance</p>
						<p class="font-semibold text-slate-800">{results.distance}</p>
					</div>
					<div>
						<p class="text-xs uppercase tracking-wide text-slate-500">Job ID</p>
						<p class="font-mono text-sm text-slate-800">{results.job_id}</p>
					</div>
				</div>

				<!-- Results Table -->
				<div class="overflow-x-auto">
					<table class="w-full text-sm">
						<thead>
							<tr class="border-b border-slate-200">
								<th class="text-left py-3 px-4 font-semibold text-slate-600">Sample</th>
								<th class="text-left py-3 px-4 font-semibold text-slate-600">Haplogroup</th>
								<th class="text-left py-3 px-4 font-semibold text-slate-600">Quality</th>
								<th class="text-right py-3 px-4 font-semibold text-slate-600">Found</th>
								<th class="text-right py-3 px-4 font-semibold text-slate-600">Missing</th>
								<th class="text-right py-3 px-4 font-semibold text-slate-600">Remaining</th>
								<th class="w-10"></th>
							</tr>
						</thead>
						<tbody>
							{#each results.results as sample (sample.id)}
								<tr class="border-b border-slate-100 hover:bg-slate-50">
									<td class="py-3 px-4 font-mono text-xs">{sample.id}</td>
									<td class="py-3 px-4 font-semibold text-blue-600">{sample.haplogroup}</td>
									<td class="py-3 px-4">
										<span class="inline-block px-2 py-1 rounded-full text-xs font-semibold text-white {getQualityClasses(sample.quality)}">
											{formatPercent(sample.quality)}
										</span>
									</td>
									<td class="py-3 px-4 text-right text-emerald-600 font-medium">{sample.top_result.found_count}</td>
									<td class="py-3 px-4 text-right text-red-500 font-medium">{sample.top_result.missing_count}</td>
									<td class="py-3 px-4 text-right text-amber-500 font-medium">{sample.top_result.remaining_count}</td>
									<td class="py-3 px-4">
										<button
											onclick={() => toggleExpanded(sample.id)}
											class="p-1 hover:bg-slate-200 rounded transition-colors"
											aria-label={expandedSamples.has(sample.id) ? 'Collapse details' : 'Expand details'}
										>
											<svg
												class="w-5 h-5 text-slate-500 transition-transform {expandedSamples.has(sample.id) ? 'rotate-180' : ''}"
												viewBox="0 0 24 24"
												fill="none"
												stroke="currentColor"
												stroke-width="2"
											>
												<polyline points="6,9 12,15 18,9" />
											</svg>
										</button>
									</td>
								</tr>
								{#if expandedSamples.has(sample.id)}
									<tr class="bg-slate-50">
										<td colspan="7" class="p-4">
											<div class="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
												<div>
													<h4 class="text-xs uppercase tracking-wide text-slate-500 mb-1">Range</h4>
													<p class="font-mono text-sm">{sample.range || 'N/A'}</p>
												</div>
												<div class="md:col-span-2 lg:col-span-3">
													<h4 class="text-xs uppercase tracking-wide text-slate-500 mb-1">Sample Polymorphisms</h4>
													<p class="font-mono text-xs text-slate-700 break-all">{sample.polymorphisms.join(' ') || 'None'}</p>
												</div>
												<div>
													<h4 class="text-xs uppercase tracking-wide text-slate-500 mb-1">Found</h4>
													<p class="font-mono text-xs text-emerald-600 break-all">{sample.top_result.found_polymorphisms.join(' ') || 'None'}</p>
												</div>
												<div>
													<h4 class="text-xs uppercase tracking-wide text-slate-500 mb-1">Missing</h4>
													<p class="font-mono text-xs text-red-500 break-all">{sample.top_result.missing_polymorphisms.join(' ') || 'None'}</p>
												</div>
												<div>
													<h4 class="text-xs uppercase tracking-wide text-slate-500 mb-1">Remaining</h4>
													<p class="font-mono text-xs text-amber-500 break-all">{sample.top_result.remaining_polymorphisms.join(' ') || 'None'}</p>
												</div>
												{#if sample.other_results.length > 0}
													<div>
														<h4 class="text-xs uppercase tracking-wide text-slate-500 mb-1">Alternative Hits</h4>
														<div class="flex flex-wrap gap-2">
															{#each sample.other_results as other, i (i)}
																<span class="inline-block px-2 py-1 bg-slate-200 rounded text-xs">
																	{other.haplogroup} ({formatPercent(other.quality)})
																</span>
															{/each}
														</div>
													</div>
												{/if}
											</div>
										</td>
									</tr>
								{/if}
							{/each}
						</tbody>
					</table>
				</div>
			</div>
		{/if}
</main>
