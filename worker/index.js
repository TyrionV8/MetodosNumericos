// Serve os arquivos do compêndio de métodos numéricos, rotinas Octave (.m)
// e recursos estáticos quando chamado diretamente ou via Service Binding pelo HardCoded.
export default {
  async fetch(request, env) {
    if (request.method !== 'GET' && request.method !== 'HEAD') {
      return new Response('Método não permitido.', { status: 405, headers: { Allow: 'GET, HEAD' } });
    }
    const response = await env.ASSETS.fetch(request);
    const headers = new Headers(response.headers);
    headers.set('Cache-Control', 'no-cache');
    headers.set('X-Content-Type-Options', 'nosniff');

    // Assegura que códigos de rotinas matemáticas sejam entregues como texto simples legível
    const url = new URL(request.url);
    if (url.pathname.endsWith('.m')) {
      headers.set('Content-Type', 'text/plain; charset=utf-8');
    }

    return new Response(response.body, { status: response.status, headers });
  },
};
