/**
 * @name shuffle
 *
 * @description
 * Fischer Yates shuffle, shuffles an existing array, maintaining the reference
 *
 * @example {{myArray | shuffle:shuffleSwitch}}
 *
 * @param filterOn {Boolean} if filterOn === false then no filtering will be performed
 *
 * @return {Array} The same array, shuffled
 *
 */
angular.module('clotho.interface').filter('shuffle',  function() {
	return function(array, filterOn) {
		if (filterOn === false || !angular.isArray(items)) {
			return items;
		} else {
			var m = array.length, t, i;

			while (m) {
				i = Math.floor(Math.random() * m--); //pick an element
				t = array[m]; //swap it with the current one
				array[m] = array[i];
				array[i] = t;
			}

			return array;
		}
	};
});