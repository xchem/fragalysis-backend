"""Project-wide Django middleware."""


class CrossOriginEmbedderPolicyMiddleware:
    """Adds 'Cross-Origin-Embedder-Policy: require-corp' to every response.

    The frontend's Moorhen viewer uses threaded WebAssembly, which the browser
    only permits when the page is cross-origin isolated. That needs COEP and
    'Cross-Origin-Opener-Policy: same-origin' on the document. Django's
    SecurityMiddleware already provides COOP (SECURE_CROSS_ORIGIN_OPENER_POLICY),
    so only COEP is added here.

    Place it first in MIDDLEWARE so it sees every response, including errors.
    """

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        response = self.get_response(request)
        response.setdefault("Cross-Origin-Embedder-Policy", "require-corp")
        return response
