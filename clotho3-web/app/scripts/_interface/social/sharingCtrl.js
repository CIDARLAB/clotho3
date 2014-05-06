angular.module('clotho.interface')
.controller('SharingCtrl', function($scope, $location, $window){

	function stdIconUrl (name) {
		return 'images/socialmedia/' + name + '.png';
	}

	$scope.social = [
		{
			"name" : "facebook",
			"prefix" : "http://www.facebook.com/sharer.php?u=",
			"icon" : stdIconUrl('facebook')
		},
		{
			"name" : "google",
			"prefix" : "https://plus.google.com/share?url=",
			"icon" : stdIconUrl('google')
		},
		{
			"name" : "twitter",
			"prefix" : "http://twitter.com/share?url=",
			"icon" : stdIconUrl('twitter')
		},
		{
			"name" : "linkedin",
			"prefix" : "http://www.linkedin.com/shareArticle?mini=true&url=",
			"icon" : stdIconUrl('linkedin')
		},
		{
			"name" : "digg",
			"prefix" : "http://www.digg.com/submit?url=",
			"icon" : stdIconUrl('digg')
		},
		{
			"name" : "reddit",
			"prefix" : "http://reddit.com/submit?url=",
			"icon" : stdIconUrl('reddit')
		},
		{
			"name" : "email",
			"prefix" : "mailto:?Body=",
			"icon" : stdIconUrl('email')
		}
	];

	$scope.share = function (site, url) {
		var shareUrl = site.prefix + ( url ? url : $location.absUrl()) ;
		$window.open(shareUrl, (site.name == 'email' ? '_self' : "_blank") );
	};

});