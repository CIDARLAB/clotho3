'use strict';

/**
 * @name unique
 * @example {{items | unique:attribute | json}}
 * Filters out all duplicate items from an array by checking the specified key
 * @param [key] {string} the name of the attribute of each object to compare for uniqueness
 if the key is empty, the entire object will be compared
 if the key === false then no filtering will be performed
 * @return {array}
 */
Application.Interface.filter('unique', ['$parse', function ($parse) {

    return function (items, filterOn) {

        if (filterOn === false) {
            return items;
        }

        if ((filterOn || angular.isUndefined(filterOn)) && angular.isArray(items)) {
            var hashCheck = {}, newItems = [],
                get = angular.isString(filterOn) ? $parse(filterOn) : function (item) { return item; };

            var extractValueToCompare = function (item) {
                return angular.isObject(item) ? get(item) : item;
            };

            angular.forEach(items, function (item) {
                var valueToCheck, isDuplicate = false;

                for (var i = 0; i < newItems.length; i++) {
                    if (angular.equals(extractValueToCompare(newItems[i]), extractValueToCompare(item))) {
                        isDuplicate = true;
                        break;
                    }
                }
                if (!isDuplicate) {
                    newItems.push(item);
                }

            });
            items = newItems;
        }
        return items;
    };
}]);

/**
 * @name highlight
 * @example

 <pre class="prettyprint">
 &lt;label&gt;&lt;input type=&quot;checkbox&quot; ng-model=&quot;caseSensitive&quot;&gt; Case Sensitive?&lt;/label&gt;
 &lt;input placeholder=&quot;Custom Class Name&quot; value=&quot;ui-match&quot; ng-model=&quot;className&quot;&gt;
 &lt;input placeholder=&quot;Enter some text to highlight&quot; value=&quot;you&quot; ng-model=&quot;highlightText&quot;&gt;
 &lt;p ng-bind-html-unsafe=&quot;&#x27;Hello there, how are you today? I\'m fine thank you.&#x27; | highlight:highlightText:className:caseSensitive&quot;&gt;&lt;/p&gt;

 &lt;style&gt;
 .ui-match { background: yellow; }
 &lt;/style&gt;
 </pre>

 * @param text {string} haystack to search through
 * @param search {string} needle to search for
 * @param {string} className Class name to use. Defaults to 'ui-match'
 * @param [caseSensitive] {boolean} optional boolean to use case-sensitive searching
 */
Application.Interface.filter('highlight', function () {
    return function (text, search, className, caseSensitive) {
        if (search || angular.isNumber(search)) {
            className = angular.isDefined(className) ? className : 'ui-match';
            text = text.toString();
            search = search.toString();
            if (caseSensitive) {
                return text.split(search).join('<span class="'+className+'">' + search + '</span>');
            } else {
                return text.replace(new RegExp(search, 'gi'), '<span class="'+className+'">$&</span>');
            }
        } else {
            return text;
        }
    };
});

/**
 * @name shuffle
 *
 * @example {{myArray | shuffle}}
 *
 * @param [key] {string} if the key === false then no filtering will be performed
 *
 * @return {array}
 *
 */
Application.Interface.filter('shuffle',  function() {
    return function (items, filterOn) {
        if (filterOn === false) { return items; }

        if ((filterOn || angular.isUndefined(filterOn)) && angular.isArray(items)) {
            var o = items.slice(0, items.length); // copy
            for(var j, x, i = o.length; i; j = parseInt(Math.random() * i), x = o[--i], o[i] = o[j], o[j] = x);
            items = o;
        }
        return items;
    }
});


/**
 * todo - major refactor
 * @name categorize
 *
 * @description
 * Organize an array by a field
 *
 * @param [key] {string} the name of the attribute of each object to compare for uniqueness
 if the key is empty, the entire object will be compared
 if the key === false then no filtering will be performed
 *
 * @returns {object} Object with keys matching param used to sort, each containing an array
 */
Application.Interface.filter('categorize', ['$parse', function($parse) {
    return function (items, filterOn) {
        if (filterOn === false) {
            return items;
        }

        if ((filterOn || angular.isUndefined(filterOn)) && angular.isArray(items)) {
            var newItems = {},
                get = angular.isString(filterOn) ? $parse(filterOn) : function (item) { return item; };

            var extractValueToCompare = function (item) {
                return angular.isObject(item) ? get(item) : item;
            };

            angular.forEach(items, function (item) {
                var type = extractValueToCompare(item);
                if (!newItems[type])
                    newItems[type] = [];
                newItems[type].push(item);
            });

            items = newItems;
        }
        return items;
    }
}]);

Application.Interface.filter('breakLines', function() {
    return function(input, charCount, joiner) {
        //todo - ignore HTML
        return (input.match(new RegExp('.{1,'+charCount+'}', 'gi')) || []).join(joiner || '\n');
    }
});

Application.Interface.filter('plainText', function() {
    // note -- not perfect, but much faster than jQuery.text()
    // http://jsperf.com/regex-replace-vs-jquery-text
    return function(input, filterOn) {
        return (filterOn) ? input.replace(/(<([^>]+)>)/ig, '') : input;
    }
});
