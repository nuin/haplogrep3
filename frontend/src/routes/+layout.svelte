<script lang="ts">
	import '../app.css';
	import favicon from '$lib/assets/favicon.svg';
	import { page } from '$app/stores';

	let { children } = $props();
</script>

<svelte:head>
	<link rel="icon" href={favicon} />
</svelte:head>

<div class="min-h-screen bg-slate-50">
	<!-- Header with Navigation -->
	<header class="bg-gradient-to-r from-blue-600 to-blue-700 text-white">
		<div class="max-w-6xl mx-auto px-4 py-6">
			<div class="flex items-center justify-between">
				<div class="flex items-center gap-3">
					<svg class="w-10 h-10" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
						<circle cx="12" cy="12" r="10" />
						<path d="M12 2a15.3 15.3 0 0 1 4 10 15.3 15.3 0 0 1-4 10 15.3 15.3 0 0 1-4-10 15.3 15.3 0 0 1 4-10z" />
						<path d="M2 12h20" />
					</svg>
					<div>
						<h1 class="text-2xl font-bold">mtclassify</h1>
						<p class="text-blue-100 text-sm">mtDNA Analysis Suite</p>
					</div>
				</div>

				<!-- Navigation -->
				<nav class="flex gap-1">
					<a
						href="/"
						class="px-4 py-2 rounded-lg text-sm font-medium transition-colors
							{$page.url.pathname === '/' ? 'bg-white/20 text-white' : 'text-blue-100 hover:bg-white/10 hover:text-white'}"
					>
						Classification
					</a>
					<a
						href="/mitomaster"
						class="px-4 py-2 rounded-lg text-sm font-medium transition-colors
							{$page.url.pathname.startsWith('/mitomaster') ? 'bg-white/20 text-white' : 'text-blue-100 hover:bg-white/10 hover:text-white'}"
					>
						MitoMaster
					</a>
				</nav>
			</div>
		</div>
	</header>

	{@render children()}
</div>
