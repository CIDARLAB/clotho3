'use strict';

angular.module('clotho.webapp').controller('BrowserCtrl',
function($scope, Clotho, $filter) {

	Clotho.recent().then(function(result) {
		$scope.recent_array = result;
		$scope.sort(false);
	});

	$scope.nav = {
		"user" : {
			"text" : "My Stuff",
			"subtext" : "Your Digital Locker",
			"value" : "user"
		},
		"group" : {
			"text" : "Our Group",
			"subtext" : "Your Group's Info",
			"value" : "group"
		},
		"all" : {
			"text" : "Everything",
			"subtext" : "All of Clotho",
			"value" : "all"
		}
	};

	$scope.setCurrent = function(value) {
		$scope.current = value;
	};

	$scope.sort = function(byCat) {
		if (byCat) {
			if (!!$scope.catSort) {return;}

			$scope.catSort = true;
			$scope.recent = $filter('categorize')($scope.recent_array, 'type');
			$scope.recent['Instance'] = $filter('categorize')($scope.recent['Instance'], 'schema.name');
		} else {
			if (!$scope.catSort && angular.isDefined($scope.catSort)) {return;}

			$scope.catSort = false;
			$scope.recent = {
				"entries" : angular.copy($scope.recent_array)
			}
		}
	};

	$scope.base64icon = base64icon;
});

//todo - move into interface service
var base64icon = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAACtElEQVR4Xu2Y3UtqURDFlyYVqQhBQon5koqYiUIElSD951qa+AkGUWDUgxC9KEhqfnfXgOLl3gKP2Xlw5kXO0T2zZ+2Z+eG2NBqNCdbYLCqAVoC2gM6ANZ6B0CGoFFAKKAWUAkqBNVZAMagYVAwqBhWDawwB/TOkGFQMKgYVg4pBxeAaK7A0Bh8fH/H6+orJZIL9/X0Eg0FYLJaZpH8og/v7e+zs7CAWi/313Ve6r8LnV7GWEuDh4QH1eh02m038D4dD+Hw++P1+eaYoqVRK3m9ubuLy8hJWq/XbeluFz+8CGhZgNBohnU7LiSYSCQwGA7y8vMDpdMLj8UjMWq0m72jb29u4uLhAu91GtVrFxsYG4vE43t/fwaS3trYQiURwc3OzsM/5ilu0mw0L0O/3ZbPj8ViS/vj4kBYIBAKyh263i0wmg8PDQ7RaLXQ6HRGKmy2VSmg2m9jd3ZV1/I5Vw/VGfS6a+PT3hgVgb5fLZfHDFmCZ05hEOByWJJl4MplEPp9Hr9ebCTAvHte4XC6cnp5iGZ+/LgBPOJvNSumen5/LKeZyOen1k5MTFItF2O12eL1ePD09SaUcHR3JM+35+Vne087OzqSKlvVpRATDFcAT5wxgb1OA6TMFCIVCqFQq/+yHA/Dq6gqcH9fX17Oq2dvbQzQanfkw4tNI8lxjWABOeFYAT+3g4EA+2ddutxvHx8dgmdOYdKFQkOR40kzu7u4Ob29vcDgcUjmsDg5ArjXq89cFYEAOsNvbWzlRGtuBvcwk542tQQFIAc4FCjKlB4Ug9zlHpjRZ1KcpFJhPkCij8UR/ylbh8397M9wCP5Wo2X5UAL0R0hshvRHSGyGzJ7GZ8ZUCSgGlgFJAKWDmFDY7tlJAKaAUUAooBcyexGbGVwooBZQCSgGlgJlT2OzYSoF1p8AnDSiNnx2jBucAAAAASUVORK5CYII=";
