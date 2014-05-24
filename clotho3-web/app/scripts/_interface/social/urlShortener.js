angular.module('clotho.urlShorten', [])
/**
 * Creates a shortened URL via google
 *
 * @example
 *
 * UrlShorten('http://www.facebook.com/myProfile#hash?query=yes').then(function (url) {
 *  console.log(url);
 * });
 */
.factory('UrlShorten',  function UrlShorten($q, $http) {
		return function (url) {
			if (angular.isUndefined(url)) {
				return $q.when(null);
			}
			return $http.post('https://www.googleapis.com/urlshortener/v1/url', {
				longUrl : url
			})
			.then(function (data) {
				return data.data.id;
			});
		}
	});